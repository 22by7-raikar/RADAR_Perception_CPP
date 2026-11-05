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
      power_(radar.samples_per_chirp, radar.chirps_per_frame) {}

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

  const std::vector<CfarCell>* candidates = &cells_;
  if (cfg_.peak_grouping) {
    suppress_non_maxima(power_, cells_, peaks_);
    candidates = &peaks_;
  }
  const auto t4 = Clock::now();

  const std::size_t nr = power_.rows();
  const std::size_t nd = power_.cols();
  for (const CfarCell& c : *candidates) {
    const std::size_t r = c.range_bin;
    const std::size_t d = c.doppler_bin;

    double dr = 0.0;
    if (r > 0 && r + 1 < nr) {
      dr = log_parabolic_offset(power_(r - 1, d), power_(r, d), power_(r + 1, d));
    }
    const double dd =
        log_parabolic_offset(power_(r, (d + nd - 1) % nd), power_(r, d), power_(r, (d + 1) % nd));

    const double range = rd_.bin_to_range(static_cast<double>(r) + dr);
    if (range < cfg_.min_range_m) continue;

    const AngleEstimate ang = angle_.estimate(spectrum_, r, d);
    Detection det;
    det.range_m = static_cast<float>(range);
    det.velocity_mps = static_cast<float>(rd_.bin_to_velocity(static_cast<double>(d) + dd));
    det.azimuth_rad = static_cast<float>(ang.azimuth_rad);
    det.x_m = static_cast<float>(range * std::cos(ang.azimuth_rad));
    det.y_m = static_cast<float>(range * std::sin(ang.azimuth_rad));
    det.snr_db =
        static_cast<float>(10.0 * std::log10(static_cast<double>(c.power) /
                                             std::max(1e-30, static_cast<double>(c.noise))));
    det.range_bin = r;
    det.doppler_bin = d;
    result.detections.push_back(det);
  }
  const auto t5 = Clock::now();

  result.timings.range_doppler_ms = elapsed_ms(t0, t1);
  result.timings.integration_ms = elapsed_ms(t1, t2);
  result.timings.cfar_ms = elapsed_ms(t2, t3);
  result.timings.peak_grouping_ms = elapsed_ms(t3, t4);
  result.timings.angle_ms = elapsed_ms(t4, t5);
  result.timings.total_ms = elapsed_ms(t0, t5);
}

}  // namespace radar
