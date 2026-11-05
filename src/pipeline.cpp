#include "radar/pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "radar/interpolation.hpp"

namespace radar {

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

CfarConfig with_looks(CfarConfig cfg, std::size_t looks) {
  cfg.num_looks = looks;
  return cfg;
}

}  // namespace

RadarPipeline::RadarPipeline(const RadarConfig& radar, const PipelineConfig& cfg)
    : radar_(radar),
      cfg_(cfg),
      rd_(radar, cfg.range_doppler),
      cfar_(with_looks(cfg.cfar, radar.num_rx)),
      angle_(radar, cfg.angle),
      spectrum_(radar.num_rx, radar.chirps_per_frame, radar.samples_per_chirp),
      power_(radar.samples_per_chirp, radar.chirps_per_frame),
      steering_(radar.num_rx) {}

void RadarPipeline::process(const DataCube& adc, FrameResult& result) {
  result.detections.clear();
  const auto t0 = Clock::now();

  rd_.process(adc, spectrum_);
  const auto t1 = Clock::now();

  integrate_power(spectrum_, power_);
  const auto t2 = Clock::now();

  cfar_.detect(power_, cells_);
  result.cfar_cells = cells_.size();
  const auto t3 = Clock::now();

  // Angle peaks for every CFAR cell. Up to max_peaks reflectors per cell.
  candidates_.clear();
  cell_first_.clear();
  for (const CfarCell& c : cells_) {
    cell_first_.push_back(candidates_.size());
    angle_.estimate_peaks(spectrum_, c.range_bin, c.doppler_bin, angles_);
    double peak_sum = 0.0;
    for (const AngleEstimate& a : angles_) peak_sum += a.peak_power;
    for (const AngleEstimate& a : angles_) {
      // Cell power is shared between angle peaks in proportion to their strength.
      candidates_.push_back({static_cast<std::size_t>(&c - cells_.data()), a.azimuth_rad,
                             static_cast<double>(c.power) * a.peak_power / peak_sum});
    }
  }
  cell_first_.push_back(candidates_.size());
  const auto t4 = Clock::now();

  const std::size_t nr = power_.rows();
  const std::size_t nd = power_.cols();
  // 3D peak grouping: a candidate survives unless a neighbouring range-Doppler
  // cell holds a stronger candidate at a similar azimuth. Plain 2D suppression
  // would discard a weak reflector whose cell borders a stronger one at a
  // different angle (a pedestrian beside a pole).
  if (cell_lookup_.rows() != nr || cell_lookup_.cols() != nd) {
    cell_lookup_.resize(nr, nd);
    cell_lookup_.fill(-1);
  }
  for (std::size_t i = 0; i < cells_.size(); ++i) {
    cell_lookup_(cells_[i].range_bin, cells_[i].doppler_bin) = static_cast<int>(i);
  }
  keep_.assign(candidates_.size(), 1);
  if (cfg_.peak_grouping) {
    for (std::size_t ci = 0; ci < candidates_.size(); ++ci) {
      const Candidate& a = candidates_[ci];
      const CfarCell& cell = cells_[a.cell];
      for (int dr = -1; dr <= 1 && keep_[ci]; ++dr) {
        const auto rr = static_cast<std::ptrdiff_t>(cell.range_bin) + dr;
        if (rr < 0 || rr >= static_cast<std::ptrdiff_t>(nr)) continue;
        for (int dd = -1; dd <= 1 && keep_[ci]; ++dd) {
          if (dr == 0 && dd == 0) continue;
          const std::size_t dc =
              static_cast<std::size_t>(static_cast<std::ptrdiff_t>(cell.doppler_bin + nd) + dd) %
              nd;
          const int other = cell_lookup_(static_cast<std::size_t>(rr), dc);
          if (other < 0) continue;
          const auto o = static_cast<std::size_t>(other);
          for (std::size_t cj = cell_first_[o]; cj < cell_first_[o + 1]; ++cj) {
            const Candidate& b = candidates_[cj];
            if (std::abs(b.azimuth_rad - a.azimuth_rad) > cfg_.grouping_angle_rad) continue;
            // Strict on one side so equal neighbours still yield one survivor.
            if (b.power > a.power || (b.power == a.power && (dr < 0 || (dr == 0 && dd < 0)))) {
              keep_[ci] = 0;
              break;
            }
          }
        }
      }
    }
  }
  for (const CfarCell& c : cells_) cell_lookup_(c.range_bin, c.doppler_bin) = -1;

  // SNR against the global floor, not the local CFAR estimate: near the
  // zero-Doppler clutter ridge the local estimate is inflated by neighbouring
  // static reflectors, which would bias RCS low.
  result.noise_floor = estimate_noise_floor(power_, radar_.num_rx, scratch_);
  const double noise = std::max(1e-30, result.noise_floor);
  for (std::size_t ci = 0; ci < candidates_.size(); ++ci) {
    if (!keep_[ci]) continue;
    const Candidate& cand = candidates_[ci];
    const std::size_t r = cells_[cand.cell].range_bin;
    const std::size_t d = cells_[cand.cell].doppler_bin;

    // Refine range and Doppler on the beam steered at this peak's azimuth,
    // not on the channel-summed power: two reflectors sharing a cell then each
    // get their own sub-bin position. The beam maximum may sit one cell away
    // from the CFAR cell when a stronger neighbour dominated the summed power.
    const double sin_az = std::sin(cand.azimuth_rad);
    const double phase_step = 2.0 * kPi * radar_.rx_spacing_wavelengths * sin_az;
    for (std::size_t k = 0; k < radar_.num_rx; ++k) {
      steering_[k] = std::polar(1.0f, static_cast<float>(phase_step * static_cast<double>(k)));
    }
    auto beam = [&](std::size_t rr, std::size_t dc) {
      cfloat acc{};
      for (std::size_t k = 0; k < radar_.num_rx; ++k) acc += spectrum_(k, dc, rr) * steering_[k];
      return static_cast<double>(squared_magnitude(acc));
    };
    auto wrap_d = [nd](std::size_t dc, int step) {
      return static_cast<std::size_t>(static_cast<std::ptrdiff_t>(dc + nd) + step) % nd;
    };
    std::size_t dp = d;
    double b0 = beam(r, dp), bl = beam(r, wrap_d(dp, -1)), bh = beam(r, wrap_d(dp, 1));
    if (bl > b0 || bh > b0) {
      dp = wrap_d(dp, bl > bh ? -1 : 1);
      b0 = beam(r, dp);
      bl = beam(r, wrap_d(dp, -1));
      bh = beam(r, wrap_d(dp, 1));
    }
    const double dd = log_parabolic_offset(bl, b0, bh);
    double dr = 0.0;
    if (r > 0 && r + 1 < nr)
      dr = log_parabolic_offset(beam(r - 1, dp), beam(r, dp), beam(r + 1, dp));
    double doppler_bin = static_cast<double>(dp) + dd;
    // Keep the velocity on the same side of the wrap as the CFAR cell.
    if (doppler_bin - static_cast<double>(d) > static_cast<double>(nd) / 2.0)
      doppler_bin -= static_cast<double>(nd);
    if (static_cast<double>(d) - doppler_bin > static_cast<double>(nd) / 2.0)
      doppler_bin += static_cast<double>(nd);
    const double range = rd_.bin_to_range(static_cast<double>(r) + dr);
    if (range < cfg_.min_range_m) continue;
    const double snr_db = 10.0 * std::log10(cand.power / noise);
    if (snr_db < cfg_.min_snr_db) continue;

    Detection det;
    det.range_m = static_cast<float>(range);
    det.velocity_mps = static_cast<float>(rd_.bin_to_velocity(doppler_bin));
    det.azimuth_rad = static_cast<float>(cand.azimuth_rad);
    det.x_m = static_cast<float>(range * std::cos(cand.azimuth_rad));
    det.y_m = static_cast<float>(range * std::sin(cand.azimuth_rad));
    det.snr_db = static_cast<float>(snr_db);
    det.rcs_dbsm = static_cast<float>(static_cast<double>(det.snr_db) + 40.0 * std::log10(range) -
                                      cfg_.rcs_offset_db);
    det.range_bin = r;
    det.doppler_bin = d;
    result.detections.push_back(det);
  }
  const auto t5 = Clock::now();

  result.timings.range_doppler_ms = elapsed_ms(t0, t1);
  result.timings.integration_ms = elapsed_ms(t1, t2);
  result.timings.cfar_ms = elapsed_ms(t2, t3);
  result.timings.angle_ms = elapsed_ms(t3, t4);
  result.timings.peak_grouping_ms = elapsed_ms(t4, t5);
  result.timings.total_ms = elapsed_ms(t0, t5);
}

double estimate_noise_floor(const Grid2D<float>& power, std::size_t num_looks,
                            std::vector<float>& scratch) {
  if (power.size() == 0) return 0.0;
  scratch.assign(power.data(), power.data() + power.size());
  auto mid = scratch.begin() + static_cast<std::ptrdiff_t>(scratch.size() / 2);
  std::nth_element(scratch.begin(), mid, scratch.end());
  const double k = static_cast<double>(std::max<std::size_t>(1, num_looks));
  const double median_over_mean = std::pow(1.0 - 1.0 / (9.0 * k), 3.0);  // Wilson-Hilferty
  return static_cast<double>(*mid) / median_over_mean;
}

}  // namespace radar
