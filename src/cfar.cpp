#include "radar/cfar.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace radar {

CaCfar2D::CaCfar2D(const CfarConfig& cfg) : cfg_(cfg) {
  if (!(cfg_.pfa > 0.0 && cfg_.pfa < 1.0))
    throw std::invalid_argument("CFAR: pfa must be in (0,1)");
  if (cfg_.num_looks == 0) throw std::invalid_argument("CFAR: num_looks must be >= 1");
  if (cfg_.train_range == 0 && cfg_.train_doppler == 0) {
    throw std::invalid_argument("CFAR: need training cells in at least one dimension");
  }
}

double CaCfar2D::false_alarm_probability(double t, std::size_t num_train, std::size_t num_looks) {
  // Cell X ~ Gamma(K), training sum S ~ Gamma(N K). For integer K:
  //   P(X > t S) = sum_{i<K} C(M+i-1, i) t^i / (1+t)^(M+i),  M = N K.
  const double m = static_cast<double>(num_train * num_looks);
  const double log_t = std::log(t);
  const double log_1pt = std::log1p(t);
  double p = 0.0;
  for (std::size_t i = 0; i < num_looks; ++i) {
    const double di = static_cast<double>(i);
    const double log_term = std::lgamma(m + di) - std::lgamma(m) - std::lgamma(di + 1.0) +
                            (i == 0 ? 0.0 : di * log_t) - (m + di) * log_1pt;
    p += std::exp(log_term);
  }
  return p;
}

double CaCfar2D::threshold_multiplier(std::size_t num_train, std::size_t num_looks, double pfa) {
  if (num_train == 0) throw std::invalid_argument("CFAR: no training cells");
  double lo = 0.0;
  double hi = 1.0;
  while (false_alarm_probability(hi, num_train, num_looks) > pfa) hi *= 2.0;
  for (int it = 0; it < 200 && (hi - lo) > 1e-12 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (false_alarm_probability(mid, num_train, num_looks) > pfa)
      lo = mid;
    else
      hi = mid;
  }
  return hi;
}

double CaCfar2D::multiplier_for(std::size_t num_train) {
  if (num_train >= multiplier_cache_.size()) {
    multiplier_cache_.resize(num_train + 1, std::numeric_limits<double>::quiet_NaN());
  }
  double& t = multiplier_cache_[num_train];
  if (std::isnan(t)) t = threshold_multiplier(num_train, cfg_.num_looks, cfg_.pfa);
  return t;
}

void CaCfar2D::detect(const Grid2D<float>& power, std::vector<CfarCell>& out) {
  out.clear();
  const std::size_t nr = power.rows();
  const std::size_t nd = power.cols();
  const std::size_t gr = cfg_.guard_range;
  const std::size_t gd = cfg_.guard_doppler;
  const std::size_t wr = cfg_.train_range + gr;
  const std::size_t wd = cfg_.train_doppler + gd;
  if (nr < 2 * wr + 1 || nd < 2 * wd + 1) {
    throw std::invalid_argument("CFAR: window larger than the power map");
  }

  // Summed-area table over a Doppler-wrapped copy: padded column c maps to
  // Doppler bin (c - wd) mod nd.
  const std::size_t pw = nd + 2 * wd;
  if (sat_.rows() != nr + 1 || sat_.cols() != pw + 1) sat_.resize(nr + 1, pw + 1);
  for (std::size_t c = 0; c <= pw; ++c) sat_(0, c) = 0.0;
  for (std::size_t r = 0; r < nr; ++r) {
    sat_(r + 1, 0) = 0.0;
    double row_sum = 0.0;
    for (std::size_t c = 0; c < pw; ++c) {
      row_sum += power(r, (c + nd - wd) % nd);
      sat_(r + 1, c + 1) = sat_(r, c + 1) + row_sum;
    }
  }
  // Inclusive box sum, rows [r0, r1], padded cols [c0, c1].
  auto box = [this](std::size_t r0, std::size_t r1, std::size_t c0, std::size_t c1) {
    return sat_(r1 + 1, c1 + 1) - sat_(r0, c1 + 1) - sat_(r1 + 1, c0) + sat_(r0, c0);
  };

  const std::size_t outer_cols = 2 * wd + 1;
  const std::size_t inner_cols = 2 * gd + 1;
  for (std::size_t r = 0; r < nr; ++r) {
    const std::size_t r0o = r >= wr ? r - wr : 0;
    const std::size_t r1o = std::min(nr - 1, r + wr);
    const std::size_t r0i = r >= gr ? r - gr : 0;
    const std::size_t r1i = std::min(nr - 1, r + gr);
    const std::size_t num_train = (r1o - r0o + 1) * outer_cols - (r1i - r0i + 1) * inner_cols;
    const double t = multiplier_for(num_train);
    const double inv_n = 1.0 / static_cast<double>(num_train);

    for (std::size_t d = 0; d < nd; ++d) {
      const std::size_t pc = d + wd;  // padded column of the CUT
      const double noise_sum = box(r0o, r1o, pc - wd, pc + wd) - box(r0i, r1i, pc - gd, pc + gd);
      const float cut = power(r, d);
      if (static_cast<double>(cut) > t * noise_sum) {
        out.push_back({r, d, cut, static_cast<float>(noise_sum * inv_n)});
      }
    }
  }
}

void suppress_non_maxima(const Grid2D<float>& power, const std::vector<CfarCell>& cells,
                         std::vector<CfarCell>& peaks) {
  peaks.clear();
  const std::size_t nr = power.rows();
  const std::size_t nd = power.cols();
  for (const CfarCell& cell : cells) {
    const float v = power(cell.range_bin, cell.doppler_bin);
    bool is_max = true;
    for (int dr = -1; dr <= 1 && is_max; ++dr) {
      const auto rr = static_cast<std::ptrdiff_t>(cell.range_bin) + dr;
      if (rr < 0 || rr >= static_cast<std::ptrdiff_t>(nr)) continue;
      for (int dd = -1; dd <= 1; ++dd) {
        if (dr == 0 && dd == 0) continue;
        const std::size_t dc =
            static_cast<std::size_t>(static_cast<std::ptrdiff_t>(cell.doppler_bin + nd) + dd) % nd;
        const float n = power(static_cast<std::size_t>(rr), dc);
        // Strict on one side so a flat-topped pair still yields one peak.
        if (n > v || (n == v && (dr < 0 || (dr == 0 && dd < 0)))) {
          is_max = false;
          break;
        }
      }
    }
    if (is_max) peaks.push_back(cell);
  }
}

}  // namespace radar
