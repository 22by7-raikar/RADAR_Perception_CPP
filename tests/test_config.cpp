#include <gtest/gtest.h>

#include <stdexcept>

#include "radar/config.hpp"

namespace radar {

TEST(RadarConfig, DerivedQuantitiesMatchDesign) {
  RadarConfig cfg;
  EXPECT_NEAR(cfg.range_resolution(), 0.9993, 1e-4);
  EXPECT_NEAR(cfg.max_range(), 255.8, 0.1);
  EXPECT_NEAR(cfg.slope(), 1.25e13, 1e9);
  EXPECT_NEAR(cfg.sample_rate(), 21.33e6, 0.01e6);
  EXPECT_NEAR(cfg.max_velocity(), 69.5, 0.1);
  EXPECT_NEAR(cfg.velocity_resolution(), 1.087, 1e-3);
  // Velocity resolution x number of chirps spans the full +/- v_max interval.
  EXPECT_NEAR(cfg.velocity_resolution() * static_cast<double>(cfg.chirps_per_frame),
              2.0 * cfg.max_velocity(), 1e-9);
  EXPECT_NO_THROW(cfg.validate());
}

TEST(RadarConfig, RejectsInvalidValues) {
  RadarConfig cfg;
  cfg.bandwidth_hz = 0.0;
  EXPECT_THROW(cfg.validate(), std::invalid_argument);

  cfg = RadarConfig{};
  cfg.chirps_per_frame = 127;
  EXPECT_THROW(cfg.validate(), std::invalid_argument);

  cfg = RadarConfig{};
  cfg.num_rx = 0;
  EXPECT_THROW(cfg.validate(), std::invalid_argument);
}

}  // namespace radar
