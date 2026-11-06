#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "radar/angle.hpp"
#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"
#include "test_helpers.hpp"

namespace radar {

namespace {

using test::kDeg;

/// Strongest detection of a single-target frame.
Detection detect_one(const Target& t, const PipelineConfig& pc = {}, std::uint32_t seed = 7) {
  const RadarConfig cfg = RadarConfig::tdm_mimo_3tx4rx();
  SimulatorConfig sc;
  sc.seed = seed;
  Simulator sim(cfg, sc);
  RadarPipeline pipe(cfg, pc);
  FrameResult res;
  pipe.process(sim.simulate_frame({t}), res);
  EXPECT_FALSE(res.detections.empty());
  if (res.detections.empty()) return {};
  return *std::max_element(
      res.detections.begin(), res.detections.end(),
      [](const Detection& a, const Detection& b) { return a.snr_db < b.snr_db; });
}

Target make_target(double range, double velocity, double az_deg, double el_deg = 0.0) {
  Target t;
  t.range_m = range;
  t.velocity_mps = velocity;
  t.azimuth_rad = az_deg * kDeg;
  t.elevation_rad = el_deg * kDeg;
  t.rcs_dbsm = 10.0;
  return t;
}

}  // namespace

TEST(Mimo, ConfigScalesVelocityLimitsWithTransmitters) {
  const RadarConfig single;
  const RadarConfig mimo = RadarConfig::tdm_mimo_3tx4rx();
  EXPECT_EQ(mimo.num_virtual(), 12u);
  EXPECT_NEAR(mimo.max_velocity(), single.max_velocity() / 3.0, 1e-12);
  EXPECT_NEAR(mimo.max_unfolded_velocity(), single.max_velocity(), 1e-12);
  EXPECT_NEAR(mimo.velocity_resolution(), single.velocity_resolution() / 3.0, 1e-12);
  EXPECT_NEAR(mimo.frame_duration(), 3.0 * single.frame_duration(), 1e-15);

  RadarConfig bad = mimo;
  bad.tx_positions.pop_back();
  EXPECT_THROW(bad.validate(), std::invalid_argument);
}

TEST(Mimo, VirtualArrayHasAzimuthRowAndElevationPairs) {
  const VirtualArray va = analyse_virtual_array(RadarConfig::tdm_mimo_3tx4rx());
  ASSERT_EQ(va.azimuth_row.size(), 8u);  // Tx0 and Tx2 rows side by side
  EXPECT_NEAR(va.azimuth_spacing, 0.5, 1e-12);
  EXPECT_EQ(va.azimuth_row.front(), 0u);  // Tx0 Rx0
  EXPECT_EQ(va.azimuth_row.back(), 11u);  // Tx2 Rx3
  ASSERT_EQ(va.elevation_pairs.size(), 4u);
  EXPECT_NEAR(va.elevation_spacing, 0.5, 1e-12);

  RadarConfig gap = RadarConfig::tdm_mimo_3tx4rx();
  gap.tx_positions[2].y = 2.5;  // leaves a hole in the azimuth row
  EXPECT_THROW(analyse_virtual_array(gap), std::invalid_argument);
}

TEST(Mimo, SimulatorMatchesDirectTdmPhase) {
  // Independent evaluation with the TDM timing written out: Tx t's chirp of
  // loop m starts at (m * Ntx + t) * Tc.
  const RadarConfig cfg = RadarConfig::tdm_mimo_3tx4rx();
  SimulatorConfig sc;
  sc.add_noise = false;
  const Target t = make_target(52.0, -31.0, 23.0, 7.0);
  Simulator sim(cfg, sc);
  const DataCube fast = sim.simulate_frame({t});

  const double c = kSpeedOfLight, lambda = cfg.wavelength(), dt = 1.0 / cfg.sample_rate();
  const double amp = sim.echo_amplitude(t);
  const double u = std::cos(t.elevation_rad) * std::sin(t.azimuth_rad);
  const double w = std::sin(t.elevation_rad);
  double max_err = 0.0;
  for (std::size_t tx = 0; tx < cfg.num_tx; ++tx) {
    for (std::size_t rx = 0; rx < cfg.num_rx; ++rx) {
      const double y = cfg.tx_positions[tx].y + 0.5 * static_cast<double>(rx);
      const double z = cfg.tx_positions[tx].z;
      for (std::size_t m = 0; m < cfg.chirps_per_frame; m += 13) {
        const double ts = static_cast<double>(m * cfg.num_tx + tx) * cfg.chirp_period();
        for (std::size_t n = 0; n < cfg.samples_per_chirp; n += 7) {
          const double tf = static_cast<double>(n) * dt;
          const double tau =
              2.0 * (t.range_m + t.velocity_mps * (ts + tf)) / c - lambda * (y * u + z * w) / c;
          double cyc =
              cfg.carrier_hz * tau + cfg.slope() * tf * tau - 0.5 * cfg.slope() * tau * tau;
          cyc -= std::floor(cyc);
          const std::complex<double> ref = std::polar(amp, 2.0 * kPi * cyc);
          const std::complex<double> got(fast(tx * cfg.num_rx + rx, m, n));
          max_err = std::max(max_err, std::abs(got - ref) / amp);
        }
      }
    }
  }
  EXPECT_LT(max_err, 1e-4);
}

struct MimoCase {
  double velocity;
  double az_deg;
  double el_deg;
};

class MimoPipeline : public ::testing::TestWithParam<MimoCase> {};

TEST_P(MimoPipeline, RecoversVelocityAzimuthAndElevation) {
  const MimoCase c = GetParam();
  const Detection d = detect_one(make_target(45.0, c.velocity, c.az_deg, c.el_deg));
  EXPECT_NEAR(d.velocity_mps, c.velocity, 0.3);
  EXPECT_NEAR(d.azimuth_rad / kDeg, c.az_deg, 1.0);
  EXPECT_NEAR(d.elevation_rad / kDeg, c.el_deg, 1.5);
  EXPECT_NEAR(d.range_m, 45.0, 0.3);
}

// The TDM Doppler span is +/-23.2 m/s; the faster cases only work with unfolding.
INSTANTIATE_TEST_SUITE_P(Targets, MimoPipeline,
                         ::testing::Values(MimoCase{0.0, 0.0, 0.0}, MimoCase{-15.0, 30.0, 0.0},
                                           MimoCase{12.0, -40.0, 5.0}, MimoCase{40.0, -20.0, 0.0},
                                           MimoCase{-60.0, 10.0, -4.0},
                                           MimoCase{-35.0, 15.0, 10.0}));

TEST(Mimo, WithoutMotionCompensationAzimuthIsBiased) {
  // A moving target's phase advances between transmitters; uncorrected, the
  // Tx2 half of the azimuth row is rotated against the Tx0 half.
  PipelineConfig off;
  off.tdm_motion_compensation = false;
  const Target t = make_target(45.0, -15.0, 30.0);
  const Detection good = detect_one(t);
  const Detection bad = detect_one(t, off);
  EXPECT_LT(std::abs(good.azimuth_rad / kDeg - 30.0), 1.0);
  EXPECT_GT(std::abs(bad.azimuth_rad / kDeg - 30.0), 3.0);
}

TEST(Mimo, WithoutUnfoldingFastTargetsAlias) {
  PipelineConfig off;
  off.velocity_unfolding = false;
  const RadarConfig cfg = RadarConfig::tdm_mimo_3tx4rx();
  const Detection d = detect_one(make_target(45.0, 40.0, -20.0), off);
  EXPECT_EQ(d.doppler_fold, 0);
  EXPECT_NEAR(d.velocity_mps, 40.0 - 2.0 * cfg.max_velocity(), 0.3);

  const Detection u = detect_one(make_target(45.0, 40.0, -20.0));
  EXPECT_EQ(u.doppler_fold, 1);
}

}  // namespace radar

namespace radar {

TEST(Pipeline, RangeDopplerCouplingIsCorrected) {
  const RadarConfig cfg;  // single transmitter: velocity is unambiguous to +/-70 m/s
  SimulatorConfig sc;
  sc.add_noise = false;
  Target t;
  t.range_m = 60.0;
  t.velocity_mps = -50.0;
  auto range_error = [&](bool correct) {
    PipelineConfig pc;
    pc.range_doppler_correction = correct;
    Simulator sim(cfg, sc);
    RadarPipeline pipe(cfg, pc);
    FrameResult res;
    pipe.process(sim.simulate_frame({t}), res);
    EXPECT_FALSE(res.detections.empty());
    double best = 1e9;
    for (const auto& d : res.detections) best = std::min(best, std::abs(d.range_m - 60.0));
    return best;
  };
  // Uncorrected: -50 m/s * (fc / S + T_frame / 2) = -0.35 m.
  EXPECT_GT(range_error(false), 0.2);
  EXPECT_LT(range_error(true), 0.1);
}

}  // namespace radar

namespace radar {

TEST(Mimo, ElevationUncertaintyMatchesMeasuredSpread) {
  const RadarConfig cfg = RadarConfig::tdm_mimo_3tx4rx();
  SimulatorConfig sc;
  sc.seed = 11;
  Simulator sim(cfg, sc);
  RadarPipeline pipe(cfg);
  FrameResult res;
  double sq_err = 0.0, sq_pred = 0.0;
  for (int i = 0; i < 150; ++i) {
    Target t;
    t.range_m = 80.0;  // ~25 dB SNR
    t.velocity_mps = -8.0;
    t.azimuth_rad = (i % 40 - 20) * test::kDeg;
    t.elevation_rad = 3.0 * test::kDeg;
    t.rcs_dbsm = 0.0;
    pipe.process(sim.simulate_frame({t}), res);
    const Detection* best = nullptr;
    for (const auto& d : res.detections)
      if (!best || d.snr_db > best->snr_db) best = &d;
    ASSERT_NE(best, nullptr);
    const double e = best->elevation_rad - t.elevation_rad;
    sq_err += e * e;
    sq_pred += static_cast<double>(best->elevation_std_rad) * best->elevation_std_rad;
  }
  const double ratio = std::sqrt(sq_err / sq_pred);
  EXPECT_GT(ratio, 0.8);
  EXPECT_LT(ratio, 1.25);
}

TEST(Mimo, SingleTransmitterHasNoElevation) {
  SimulatorConfig sc;
  sc.add_noise = false;
  const RadarConfig cfg;
  Simulator sim(cfg, sc);
  RadarPipeline pipe(cfg);
  FrameResult res;
  pipe.process(sim.simulate_frame({Target{40.0, 5.0, 0.2, 10.0}}), res);
  ASSERT_FALSE(res.detections.empty());
  EXPECT_EQ(res.detections.front().elevation_rad, 0.0f);
  EXPECT_TRUE(std::isinf(res.detections.front().elevation_std_rad));
}

}  // namespace radar
