#include <gtest/gtest.h>

#include "radar/angle.hpp"
#include "radar/range_doppler.hpp"
#include "radar/simulator.hpp"
#include "test_helpers.hpp"

namespace radar {

class AngleTest : public ::testing::TestWithParam<double> {};

TEST_P(AngleTest, RecoversAzimuthOfNoiselessTarget) {
  const double az_deg = GetParam();
  RadarConfig cfg;
  SimulatorConfig sc;
  sc.add_noise = false;
  Simulator sim(cfg, sc);
  Target t;
  t.range_m = 40.0;
  t.velocity_mps = 3.0;
  t.azimuth_rad = az_deg * test::kDeg;

  RangeDopplerProcessor rd(cfg);
  DataCube spec;
  rd.process(sim.simulate_frame({t}), spec);
  Grid2D<float> p;
  integrate_power(spec, p);
  const auto [r, d] = test::argmax(p);

  AngleEstimator est(cfg);
  const AngleEstimate a = est.estimate(spec, r, d);
  EXPECT_NEAR(a.azimuth_rad / test::kDeg, az_deg, 1.0);
}

INSTANTIATE_TEST_SUITE_P(Azimuths, AngleTest,
                         ::testing::Values(-50.0, -30.0, -10.0, 0.0, 5.0, 20.0, 45.0));

TEST(Angle, RejectsTooSmallFft) {
  RadarConfig cfg;
  AngleConfig ac;
  ac.fft_size = cfg.num_rx - 1;
  EXPECT_THROW(AngleEstimator(cfg, ac), std::invalid_argument);
}

}  // namespace radar
