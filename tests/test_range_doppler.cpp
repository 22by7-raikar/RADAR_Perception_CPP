#include <gtest/gtest.h>

#include "radar/range_doppler.hpp"
#include "radar/simulator.hpp"
#include "test_helpers.hpp"

namespace radar {

namespace {

Grid2D<float> power_for(const RadarConfig& cfg, const Target& t) {
  SimulatorConfig sc;
  sc.add_noise = false;
  Simulator sim(cfg, sc);
  const DataCube adc = sim.simulate_frame({t});
  RangeDopplerProcessor rd(cfg);
  DataCube spec;
  rd.process(adc, spec);
  Grid2D<float> p;
  integrate_power(spec, p);
  return p;
}

}  // namespace

TEST(RangeDoppler, PeakLandsInExpectedBins) {
  RadarConfig cfg;
  const std::size_t r_bin = 50;
  const int v_bin = 10;
  Target t;
  t.range_m = static_cast<double>(r_bin) * cfg.range_resolution();
  t.velocity_mps = v_bin * cfg.velocity_resolution();

  const auto p = power_for(cfg, t);
  ASSERT_EQ(p.rows(), cfg.samples_per_chirp);
  ASSERT_EQ(p.cols(), cfg.chirps_per_frame);
  const auto [r, d] = test::argmax(p);
  EXPECT_EQ(r, r_bin);
  EXPECT_EQ(static_cast<int>(d) - static_cast<int>(cfg.chirps_per_frame / 2), v_bin);
}

TEST(RangeDoppler, ApproachingTargetHasNegativeVelocityBin) {
  RadarConfig cfg;
  Target t;
  t.range_m = 80.0;
  t.velocity_mps = -20.0;
  const auto p = power_for(cfg, t);
  const auto [r, d] = test::argmax(p);
  RangeDopplerProcessor rd(cfg);
  EXPECT_NEAR(rd.bin_to_range(static_cast<double>(r)), 80.0, cfg.range_resolution());
  EXPECT_NEAR(rd.bin_to_velocity(static_cast<double>(d)), -20.0, cfg.velocity_resolution());
}

TEST(RangeDoppler, RejectsMismatchedFrame) {
  RadarConfig cfg;
  RangeDopplerProcessor rd(cfg);
  DataCube wrong(cfg.num_rx, cfg.chirps_per_frame, cfg.samples_per_chirp / 2);
  DataCube spec;
  EXPECT_THROW(rd.process(wrong, spec), std::invalid_argument);
}

}  // namespace radar
