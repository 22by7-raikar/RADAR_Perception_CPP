#include "radar/impairments.hpp"

#include <cmath>
#include <random>

namespace radar {

void add_interference(const RadarConfig& radar, const Interferer& in, DataCube& frame) {
  const std::size_t nch = radar.num_virtual();
  const std::size_t nd = radar.chirps_per_frame;
  const std::size_t nr = radar.samples_per_chirp;
  const double fs = radar.sample_rate();
  const double dt = 1.0 / fs;
  const double amp = std::sqrt(radar.noise_power * std::pow(10.0, in.inr_db / 10.0));
  const double u = std::sin(in.azimuth_rad);

  for (std::size_t ch = 0; ch < nch; ++ch) {
    // Only the receive element matters for a one-way path.
    const double rx_phase = -2.0 * kPi * radar.rx_position(ch % radar.num_rx).y * u;
    const double tx_offset = static_cast<double>(ch / radar.num_rx) * radar.chirp_period();
    for (std::size_t m = 0; m < nd; ++m) {
      const double t0 = static_cast<double>(m) * radar.loop_period() + tx_offset;
      cfloat* row = &frame(ch, m, 0);
      // The IF phase is the running integral of the beat frequency; only
      // samples whose beat lies inside the IF band are added.
      double cycles = 0.0;
      for (std::size_t n = 0; n < nr; ++n) {
        const double t = t0 + static_cast<double>(n) * dt;
        const double since = t - in.time_offset_s;
        const double k = std::floor(since / in.period_s);
        const double in_ramp = since - k * in.period_s;
        if (since < 0.0 || in_ramp > in.ramp_duration_s) continue;  // aggressor idle
        const double f_aggressor = in.start_frequency_hz + in.slope_hz_per_s * in_ramp;
        const double f_victim = radar.carrier_hz + radar.slope() * (t - t0);
        const double beat = f_victim - f_aggressor;
        cycles += beat * dt;
        if (beat < 0.0 || beat >= fs) continue;  // outside the IF filter
        const double phase = 2.0 * kPi * (cycles - std::floor(cycles)) + rx_phase;
        row[n] += cfloat(static_cast<float>(amp * std::cos(phase)),
                         static_cast<float>(amp * std::sin(phase)));
      }
    }
  }
}

std::vector<Target> clutter_targets(const ClutterPatch& patch, std::uint32_t seed) {
  std::mt19937 rng(seed);
  // Uniform over the area, not over range: density grows with range.
  std::uniform_real_distribution<double> r2(patch.range_min_m * patch.range_min_m,
                                            patch.range_max_m * patch.range_max_m);
  std::uniform_real_distribution<double> az(patch.azimuth_min_rad, patch.azimuth_max_rad);
  std::normal_distribution<double> v(0.0, patch.velocity_spread_mps);
  std::vector<Target> out;
  out.reserve(patch.count);
  for (std::size_t i = 0; i < patch.count; ++i) {
    Target t;
    t.range_m = std::sqrt(r2(rng));
    t.azimuth_rad = az(rng);
    t.velocity_mps = v(rng);
    t.rcs_dbsm = patch.rcs_dbsm;
    out.push_back(t);
  }
  return out;
}

namespace {

double radial_velocity(const Vec2& p, const Vec2& v) { return p.dot(v) / p.norm(); }

}  // namespace

Target direct_target(const MovingPoint& p) {
  Target t;
  t.range_m = p.position.norm();
  t.azimuth_rad = std::atan2(p.position.y(), p.position.x());
  t.velocity_mps = radial_velocity(p.position, p.velocity);
  t.rcs_dbsm = p.rcs_dbsm;
  return t;
}

std::vector<Target> multipath_targets(const MovingPoint& p, const PlanarReflector& wall) {
  // Sensor and point must be on the same side of the reflector.
  if ((p.position.y() - wall.y_m) * (0.0 - wall.y_m) <= 0.0) return {};
  const Vec2 img(p.position.x(), 2.0 * wall.y_m - p.position.y());
  const Vec2 img_v(p.velocity.x(), -p.velocity.y());
  const double r_direct = p.position.norm(), r_image = img.norm();
  const double v_direct = radial_velocity(p.position, p.velocity);
  const double v_image = radial_velocity(img, img_v);

  std::vector<Target> out;
  Target dbl;
  dbl.range_m = r_image;
  dbl.azimuth_rad = std::atan2(img.y(), img.x());
  dbl.velocity_mps = v_image;
  dbl.rcs_dbsm = p.rcs_dbsm - 2.0 * wall.loss_db;
  out.push_back(dbl);

  for (const Vec2* arrival : {&p.position, &img}) {
    Target s;
    s.range_m = 0.5 * (r_direct + r_image);
    s.azimuth_rad = std::atan2(arrival->y(), arrival->x());
    s.velocity_mps = 0.5 * (v_direct + v_image);
    s.rcs_dbsm = p.rcs_dbsm - wall.loss_db;
    out.push_back(s);
  }
  return out;
}

}  // namespace radar
