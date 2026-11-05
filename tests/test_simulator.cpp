#include <gtest/gtest.h>

#include <cmath>

#include "radar/simulator.hpp"

namespace radar {

namespace {

// Direct evaluation of the dechirped phase at every sample. Slow, but it is the
// definition the fast recurrence in Simulator must reproduce.
DataCube reference_frame(const RadarConfig& cfg, const SimulatorConfig& sc,
                         const std::vector<Target>& targets) {
  Simulator amp_source(cfg, sc);
  DataCube out(cfg.num_rx, cfg.chirps_per_frame, cfg.samples_per_chirp);
  const double c = kSpeedOfLight;
  const double dt = 1.0 / cfg.sample_rate();
  const double d = cfg.rx_spacing_wavelengths * cfg.wavelength();
  for (const Target& t : targets) {
    const double amp = amp_source.echo_amplitude(t);
    for (std::size_t k = 0; k < cfg.num_rx; ++k) {
      for (std::size_t m = 0; m < cfg.chirps_per_frame; ++m) {
        for (std::size_t n = 0; n < cfg.samples_per_chirp; ++n) {
          const double tf = static_cast<double>(n) * dt;
          const double ts = static_cast<double>(m) * cfg.chirp_period();
          const double tau = 2.0 * (t.range_m + t.velocity_mps * (ts + tf)) / c -
                             static_cast<double>(k) * d * std::sin(t.azimuth_rad) / c;
          double cyc =
              cfg.carrier_hz * tau + cfg.slope() * tf * tau - 0.5 * cfg.slope() * tau * tau;
          cyc -= std::floor(cyc);
          out(k, m, n) += std::polar(static_cast<float>(amp), static_cast<float>(2.0 * kPi * cyc));
        }
      }
    }
  }
  return out;
}

}  // namespace

TEST(Simulator, FastRecurrenceMatchesDirectPhaseEvaluation) {
  RadarConfig cfg;
  SimulatorConfig sc;
  sc.add_noise = false;
  const std::vector<Target> targets{{37.3, -21.7, 0.31, 12.0}, {201.0, 55.0, -0.9, 20.0}};

  Simulator sim(cfg, sc);
  const DataCube fast = sim.simulate_frame(targets);
  const DataCube ref = reference_frame(cfg, sc, targets);

  double max_err = 0.0, max_mag = 0.0;
  for (std::size_t i = 0; i < ref.size(); ++i) {
    max_err = std::max(max_err, static_cast<double>(std::abs(fast.data()[i] - ref.data()[i])));
    max_mag = std::max(max_mag, static_cast<double>(std::abs(ref.data()[i])));
  }
  EXPECT_LT(max_err / max_mag, 1e-4);
}

TEST(Simulator, EchoAmplitudeFollowsRadarEquation) {
  RadarConfig cfg;
  Simulator sim(cfg);
  const double a100 = sim.echo_amplitude({100.0, 0.0, 0.0, 10.0});
  const double a200 = sim.echo_amplitude({200.0, 0.0, 0.0, 10.0});
  const double a100_big = sim.echo_amplitude({100.0, 0.0, 0.0, 20.0});
  EXPECT_NEAR(20.0 * std::log10(a100 / a200), 40.0 * std::log10(2.0), 1e-9);  // R^4 in power
  EXPECT_NEAR(20.0 * std::log10(a100_big / a100), 10.0, 1e-9);
  EXPECT_THROW(sim.echo_amplitude({0.0, 0.0, 0.0, 0.0}), std::invalid_argument);
}

}  // namespace radar
