#include "radar/mitigation.hpp"

#include <pocketfft_hdronly.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "radar/window.hpp"

namespace radar {

// --- Interference ------------------------------------------------------------

namespace {

float median_of(const float* data, std::size_t n, std::vector<float>& scratch) {
  scratch.assign(data, data + n);
  auto mid = scratch.begin() + static_cast<std::ptrdiff_t>(n / 2);
  std::nth_element(scratch.begin(), mid, scratch.end());
  return *mid;
}

}  // namespace

std::vector<float> notched_residual(const cfloat* samples, std::size_t n, double notch_db) {
  const std::vector<float> window = make_window(WindowType::kHann, n);
  const pocketfft::detail::pocketfft_c<float> plan(n);
  std::vector<cfloat> buf(n);
  std::vector<float> pw(n), scratch;
  for (std::size_t i = 0; i < n; ++i) buf[i] = samples[i] * window[i];
  auto* data = reinterpret_cast<pocketfft::detail::cmplx<float>*>(buf.data());
  plan.exec(data, 1.0f, true);
  for (std::size_t k = 0; k < n; ++k) pw[k] = squared_magnitude(buf[k]);
  const float notch =
      static_cast<float>(std::pow(10.0, notch_db / 10.0)) * median_of(pw.data(), n, scratch);
  for (std::size_t k = 0; k < n; ++k)
    if (pw[k] > notch) buf[k] = cfloat{};
  plan.exec(data, 1.0f / static_cast<float>(n), false);
  for (std::size_t i = 0; i < n; ++i) pw[i] = squared_magnitude(buf[i]);
  return pw;
}

std::size_t mitigate_interference(const RadarConfig& radar, const InterferenceMitigationConfig& cfg,
                                  DataCube& adc, std::vector<float>& scratch) {
  const std::size_t nrx = radar.num_rx;
  const std::size_t nd = adc.num_chirps();
  const std::size_t nr = adc.num_samples();
  const std::size_t groups = adc.num_rx() / nrx;  // one group per transmitter (time slot)
  const float notch = static_cast<float>(std::pow(10.0, cfg.notch_db / 10.0));
  const float factor = static_cast<float>(std::pow(10.0, cfg.threshold_db / 10.0));
  const auto guard = static_cast<std::ptrdiff_t>(cfg.guard_samples);
  const float inv_n = 1.0f / static_cast<float>(nr);
  // One plan, executed in place: the convenience c2c() sets up shapes and
  // buffers per call, which dominated this loop (thousands of short FFTs).
  const pocketfft::detail::pocketfft_c<float> plan(nr);
  auto fft = [&plan](std::vector<cfloat>& v, bool forward, float scale) {
    plan.exec(reinterpret_cast<pocketfft::detail::cmplx<float>*>(v.data()), scale, forward);
  };

  std::vector<cfloat> buf(nr);
  std::vector<float> pw(nr);
  std::vector<char> flag(nr);
  const std::vector<float> window = make_window(WindowType::kHann, nr);
  std::vector<float> w2(nr), scratch_center(nr);
  for (std::size_t n = 0; n < nr; ++n) w2[n] = window[n] * window[n];
  std::size_t blanked = 0;
  for (std::size_t g = 0; g < groups; ++g) {
    for (std::size_t m = 0; m < nd; ++m) {
      // Residual of the first receiver with narrowband components removed.
      // Hann-windowed so each echo stays in its main lobe and the notch takes
      // all of it (rectangular leakage from off-bin tones would remain).
      const cfloat* row = &adc(g * nrx, m, 0);
      for (std::size_t n = 0; n < nr; ++n) buf[n] = row[n] * window[n];
      fft(buf, true, 1.0f);
      for (std::size_t k = 0; k < nr; ++k) pw[k] = squared_magnitude(buf[k]);
      const float spec_median = median_of(pw.data(), nr, scratch);
      if (!(spec_median > 0.0f)) continue;
      for (std::size_t k = 0; k < nr; ++k)
        if (pw[k] > notch * spec_median) buf[k] = cfloat{};
      fft(buf, false, inv_n);
      // Undo the window's taper in the comparison: threshold each sample
      // against the residual level expected at its window weight. The weight
      // is floored at 0.05 (window 0.22) so leftover leakage near the chirp
      // edges is not amplified into false flags; a 30 dB burst there still
      // clears it. Samples where the window is below 0.1 are not tested; the
      // range window attenuates a burst there by > 20 dB anyway.
      std::size_t central = 0;
      for (std::size_t n = 0; n < nr; ++n) {
        pw[n] = squared_magnitude(buf[n]);
        if (w2[n] >= 0.25f) scratch_center[central++] = pw[n] / w2[n];
      }
      const float thr = factor * median_of(scratch_center.data(), central, scratch);

      std::fill(flag.begin(), flag.end(), 0);
      bool any = false;
      for (std::size_t n = 0; n < nr; ++n) {
        if (w2[n] < 0.01f || pw[n] <= thr * std::max(w2[n], 0.05f)) continue;
        any = true;
        const auto lo = std::max<std::ptrdiff_t>(0, static_cast<std::ptrdiff_t>(n) - guard);
        const auto hi = std::min<std::ptrdiff_t>(static_cast<std::ptrdiff_t>(nr) - 1,
                                                 static_cast<std::ptrdiff_t>(n) + guard);
        for (auto k = lo; k <= hi; ++k) flag[static_cast<std::size_t>(k)] = 1;
      }
      if (!any) continue;
      std::size_t gap = 0;
      for (std::size_t n = 0; n < nr; ++n) gap += flag[n] ? 1 : 0;
      blanked += gap;
      if (cfg.reconstruct_iterations == 0 || gap >= nr / 2) {
        for (std::size_t n = 0; n < nr; ++n) {
          if (!flag[n]) continue;
          for (std::size_t r = 0; r < nrx; ++r) adc(g * nrx + r, m, n) = cfloat{};
        }
        continue;
      }
      // Refill the gap per receiver by iterative thresholded reconstruction
      // (IMAT): keep the spectral bins above a threshold that drops each
      // iteration, transform back, and overwrite only the gap. Converges for
      // spectra dominated by a few tones, including their sidelobes.
      for (std::size_t r = 0; r < nrx; ++r) {
        cfloat* x = &adc(g * nrx + r, m, 0);
        for (std::size_t n = 0; n < nr; ++n)
          if (flag[n]) x[n] = cfloat{};
        float floor_power = 0.0f, peak_power = 0.0f;
        const std::size_t iterations = cfg.reconstruct_iterations;
        for (std::size_t it = 0; it < iterations; ++it) {
          std::copy(x, x + nr, buf.begin());
          fft(buf, true, 1.0f);
          float peak = 0.0f;
          for (std::size_t k = 0; k < nr; ++k) {
            pw[k] = squared_magnitude(buf[k]);
            peak = std::max(peak, pw[k]);
          }
          if (it == 0) {
            floor_power = notch * median_of(pw.data(), nr, scratch);
            peak_power = peak;
          }
          // Thresholds evenly spaced in dB from the first peak down to the
          // notch level, reaching it on the last iteration.
          const double frac = static_cast<double>(it + 1) / static_cast<double>(iterations);
          const float keep = std::max(
              floor_power,
              static_cast<float>(static_cast<double>(peak_power) *
                                 std::pow(static_cast<double>(floor_power) /
                                              std::max(1e-30, static_cast<double>(peak_power)),
                                          frac)));
          for (std::size_t k = 0; k < nr; ++k)
            if (!(pw[k] > keep)) buf[k] = cfloat{};
          fft(buf, false, inv_n);
          for (std::size_t n = 0; n < nr; ++n)
            if (flag[n]) x[n] = buf[n];
        }
      }
    }
  }
  return blanked;
}

// --- Ordered-statistic CFAR --------------------------------------------------

namespace {

/// P(Gamma(K, 1) > x) for integer K.
double gamma_survival(double x, std::size_t k) {
  if (x <= 0.0) return 1.0;
  double term = 1.0, sum = 1.0;
  for (std::size_t i = 1; i < k; ++i) {
    term *= x / static_cast<double>(i);
    sum += term;
  }
  return std::exp(-x) * sum;
}

double gamma_pdf(double x, std::size_t k) {
  if (x <= 0.0) return 0.0;
  const auto kd = static_cast<double>(k);
  return std::exp((kd - 1.0) * std::log(x) - x - std::lgamma(kd));
}

}  // namespace

OsCfar2D::OsCfar2D(const CfarConfig& cfg, double rank) : cfg_(cfg), rank_(rank) {
  if (cfg_.num_looks == 0) throw std::invalid_argument("OS-CFAR: num_looks must be >= 1");
  if (!(rank_ > 0.0 && rank_ <= 1.0)) throw std::invalid_argument("OS-CFAR: rank in (0, 1]");
  if (!(cfg_.pfa > 0.0 && cfg_.pfa < 1.0)) throw std::invalid_argument("OS-CFAR: pfa in (0, 1)");
}

std::size_t OsCfar2D::rank_for(std::size_t num_train) const {
  const auto k = static_cast<std::size_t>(std::lround(rank_ * static_cast<double>(num_train)));
  return std::clamp<std::size_t>(k, 1, num_train);
}

double OsCfar2D::false_alarm_probability(double t, std::size_t n, std::size_t k,
                                         std::size_t looks) {
  // Pfa = E_y[ P(X > t y) ], y the k-th smallest of n Gamma(K) samples:
  //   f_k(y) = n!/((k-1)!(n-k)!) F^(k-1) (1-F)^(n-k) f(y).
  const auto nd = static_cast<double>(n), kd = static_cast<double>(k);
  const double log_c = std::lgamma(nd + 1.0) - std::lgamma(kd) - std::lgamma(nd - kd + 1.0);
  const double kl = static_cast<double>(looks);
  const double y_max = kl + 40.0 * std::sqrt(kl) + 40.0;
  constexpr int kSteps = 8000;
  const double h = y_max / kSteps;
  double sum = 0.0;
  for (int i = 1; i < kSteps; ++i) {  // integrand is 0 at both ends
    const double y = h * i;
    const double surv = gamma_survival(y, looks);
    const double cdf = 1.0 - surv;
    if (cdf <= 0.0 || surv <= 0.0) continue;
    const double log_f = log_c + (kd - 1.0) * std::log(cdf) + (nd - kd) * std::log(surv) +
                         std::log(gamma_pdf(y, looks));
    sum += std::exp(log_f) * gamma_survival(t * y, looks);
  }
  return sum * h;
}

double OsCfar2D::threshold_multiplier(std::size_t n, std::size_t k, std::size_t looks, double pfa) {
  double lo = 0.0, hi = 1.0;
  while (false_alarm_probability(hi, n, k, looks) > pfa) hi *= 2.0;
  for (int i = 0; i < 60; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (false_alarm_probability(mid, n, k, looks) > pfa)
      lo = mid;
    else
      hi = mid;
  }
  return hi;
}

constexpr double OsCfar2D::kPrescreen;

double OsCfar2D::multiplier_for(std::size_t n) {
  if (multiplier_cache_.size() <= n) multiplier_cache_.resize(n + 1, std::nan(""));
  double& t = multiplier_cache_[n];
  if (std::isnan(t)) t = threshold_multiplier(n, rank_for(n), cfg_.num_looks, cfg_.pfa);
  return t;
}

void OsCfar2D::prepare(std::size_t range_bins) {
  const auto nr = static_cast<std::ptrdiff_t>(range_bins);
  const auto w = static_cast<std::ptrdiff_t>(cfg_.train_range + cfg_.guard_range);
  const auto gr = static_cast<std::ptrdiff_t>(cfg_.guard_range);
  const auto width = static_cast<std::ptrdiff_t>(2 * (cfg_.train_doppler + cfg_.guard_doppler) + 1);
  const auto guard_width = static_cast<std::ptrdiff_t>(2 * cfg_.guard_doppler + 1);
  auto span = [nr](std::ptrdiff_t lo, std::ptrdiff_t hi) {
    return std::max<std::ptrdiff_t>(0, std::min(hi, nr - 1) - std::max<std::ptrdiff_t>(lo, 0) + 1);
  };
  for (std::ptrdiff_t r = 0; r < nr; ++r) {
    if (r > w && r < nr - 1 - w) continue;  // interior rows share the full size
    const std::ptrdiff_t n = span(r - w, r + w) * width - span(r - gr, r + gr) * guard_width;
    if (n > 0) multiplier_for(static_cast<std::size_t>(n));
  }
  if (nr > 2 * w + 1) {
    multiplier_for(static_cast<std::size_t>((2 * w + 1) * width - (2 * gr + 1) * guard_width));
  }
}

void OsCfar2D::detect(const Grid2D<float>& power, std::vector<CfarCell>& out) {
  out.clear();
  const auto nr = static_cast<std::ptrdiff_t>(power.rows());
  const auto nd = static_cast<std::ptrdiff_t>(power.cols());
  if (nr == 0 || nd == 0) return;
  const auto tr = static_cast<std::ptrdiff_t>(cfg_.train_range);
  const auto td = static_cast<std::ptrdiff_t>(cfg_.train_doppler);
  const auto gr = static_cast<std::ptrdiff_t>(cfg_.guard_range);
  const auto gd = static_cast<std::ptrdiff_t>(cfg_.guard_doppler);
  const std::size_t full = static_cast<std::size_t>((2 * (tr + gr) + 1) * (2 * (td + gd) + 1) -
                                                    (2 * gr + 1) * (2 * gd + 1));

  // Pre-screen against the map median (see class comment). In homogeneous
  // noise the order statistic sits above the median (rank > 0.5) with a
  // spread of a few percent, so a cell below 0.8 x threshold x median could
  // only pass if the order statistic fell ~8 sigma low.
  scratch_.assign(power.data(), power.data() + power.size());
  auto mid = scratch_.begin() + static_cast<std::ptrdiff_t>(scratch_.size() / 2);
  std::nth_element(scratch_.begin(), mid, scratch_.end());
  const double pre = kPrescreen * multiplier_for(full) * static_cast<double>(*mid);

  train_.reserve(full);
  for (std::ptrdiff_t r = 0; r < nr; ++r) {
    for (std::ptrdiff_t d = 0; d < nd; ++d) {
      const float p = power(static_cast<std::size_t>(r), static_cast<std::size_t>(d));
      if (static_cast<double>(p) <= pre) continue;
      train_.clear();
      for (std::ptrdiff_t dr = -(tr + gr); dr <= tr + gr; ++dr) {
        const std::ptrdiff_t rr = r + dr;
        if (rr < 0 || rr >= nr) continue;  // clipped at the range edges
        for (std::ptrdiff_t dd = -(td + gd); dd <= td + gd; ++dd) {
          if (std::abs(dr) <= gr && std::abs(dd) <= gd) continue;  // guard + cell under test
          const std::ptrdiff_t dc = ((d + dd) % nd + nd) % nd;     // Doppler wraps
          train_.push_back(power(static_cast<std::size_t>(rr), static_cast<std::size_t>(dc)));
        }
      }
      const std::size_t n = train_.size();
      if (n == 0) continue;
      const std::size_t k = rank_for(n);
      std::nth_element(train_.begin(), train_.begin() + static_cast<std::ptrdiff_t>(k - 1),
                       train_.end());
      const float y = train_[k - 1];
      if (static_cast<double>(p) > multiplier_for(n) * static_cast<double>(y)) {
        out.push_back({static_cast<std::size_t>(r), static_cast<std::size_t>(d), p, y});
      }
    }
  }
}

// --- Multipath ghosts --------------------------------------------------------

namespace {

bool is_static(const Detection& d, const Vec2& ego_velocity, double tol) {
  const double expected =
      -(ego_velocity.x() * std::cos(d.azimuth_rad) + ego_velocity.y() * std::sin(d.azimuth_rad));
  return std::abs(static_cast<double>(d.velocity_mps) - expected) < tol;
}

}  // namespace

double ReflectorLine::sensor_side_distance(double x, double y) const {
  const double d = (y - slope * x - offset) / std::sqrt(1.0 + slope * slope);
  return offset > 0.0 ? -d : d;  // the sensor (origin) is at raw distance -offset
}

ReflectorLine fit_reflector(const std::vector<Detection>& dets, const Vec2& ego_velocity, int side,
                            const GhostFilterConfig& cfg) {
  std::vector<Vec2> pts;
  for (const auto& d : dets) {
    if (static_cast<double>(d.y_m) * side < cfg.min_lateral_m) continue;
    if (!is_static(d, ego_velocity, cfg.static_velocity_tol_mps)) continue;
    pts.emplace_back(d.x_m, d.y_m);
  }
  ReflectorLine best;
  if (pts.size() < std::max<std::size_t>(2, cfg.min_inliers)) return best;

  std::mt19937 rng(cfg.seed);
  std::uniform_int_distribution<std::size_t> pick(0, pts.size() - 1);
  auto count = [&](double s, double o) {
    const double norm = std::sqrt(1.0 + s * s);
    std::size_t c = 0;
    for (const Vec2& p : pts)
      if (std::abs(p.y() - s * p.x() - o) / norm < cfg.inlier_distance_m) ++c;
    return c;
  };
  double best_s = 0.0, best_o = 0.0;
  std::size_t best_n = 0;
  for (std::size_t it = 0; it < cfg.ransac_iterations; ++it) {
    const Vec2& a = pts[pick(rng)];
    const Vec2& b = pts[pick(rng)];
    if (std::abs(a.x() - b.x()) < 1.0) continue;
    const double s = (b.y() - a.y()) / (b.x() - a.x());
    if (std::abs(s) > cfg.max_slope) continue;
    const double o = a.y() - s * a.x();
    const std::size_t n = count(s, o);
    if (n > best_n) {
      best_n = n;
      best_s = s;
      best_o = o;
    }
  }
  if (best_n < cfg.min_inliers) return best;

  // Least-squares refit on the inliers.
  const double norm = std::sqrt(1.0 + best_s * best_s);
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  std::size_t n = 0;
  for (const Vec2& p : pts) {
    if (std::abs(p.y() - best_s * p.x() - best_o) / norm >= cfg.inlier_distance_m) continue;
    sx += p.x();
    sy += p.y();
    sxx += p.x() * p.x();
    sxy += p.x() * p.y();
    ++n;
  }
  const double nn = static_cast<double>(n);
  const double den = nn * sxx - sx * sx;
  if (std::abs(den) > 1e-9) {
    best.slope = (nn * sxy - sx * sy) / den;
    best.offset = (sy - best.slope * sx) / nn;
  } else {
    best.slope = best_s;
    best.offset = best_o;
  }
  best.inliers = n;
  best.valid = true;
  return best;
}

void flag_multipath_ghosts(const std::vector<Detection>& dets,
                           const std::vector<ReflectorLine>& reflectors, const Vec2& ego_velocity,
                           const GhostFilterConfig& cfg, std::vector<char>& ghost) {
  ghost.assign(dets.size(), 0);
  for (std::size_t i = 0; i < dets.size(); ++i) {
    const Detection& d = dets[i];
    if (is_static(d, ego_velocity, cfg.static_velocity_tol_mps)) continue;
    for (const ReflectorLine& line : reflectors) {
      if (!line.valid) continue;
      const double dist = line.sensor_side_distance(d.x_m, d.y_m);
      if (dist > -cfg.behind_margin_m) continue;
      if (!cfg.require_mirror_match) {
        ghost[i] = 1;
        break;
      }
      // Mirror the detection across the line and look for the real object.
      const double norm = std::sqrt(1.0 + line.slope * line.slope);
      const double raw = (d.y_m - line.slope * d.x_m - line.offset) / norm;
      const Vec2 normal(-line.slope / norm, 1.0 / norm);
      const Vec2 mirror = Vec2(d.x_m, d.y_m) - 2.0 * raw * normal;
      for (std::size_t j = 0; j < dets.size() && !ghost[i]; ++j) {
        if (j == i) continue;
        const Detection& o = dets[j];
        if (line.sensor_side_distance(o.x_m, o.y_m) < 0.0) continue;
        if ((Vec2(o.x_m, o.y_m) - mirror).norm() > cfg.mirror_gate_m) continue;
        if (std::abs(o.velocity_mps - d.velocity_mps) > cfg.mirror_velocity_gate_mps) continue;
        ghost[i] = 1;
      }
      if (ghost[i]) break;
    }
  }
}

}  // namespace radar
