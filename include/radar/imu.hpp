#pragma once

#include <cstdint>
#include <random>

#include "radar/geometry.hpp"
#include "radar/scene.hpp"

namespace radar {

/// Planar IMU sample: specific force along body x/y and yaw rate.
struct ImuSample {
  double time_s = 0.0;
  Vec2 accel_body = Vec2::Zero();
  double gyro_z = 0.0;
};

/// Consumer-grade MEMS IMU: white noise plus constant bias with a slow random walk.
struct ImuConfig {
  double rate_hz = 100.0;
  double accel_noise_std = 0.05;  ///< m/s^2 per sample
  double gyro_noise_std = 0.002;  ///< rad/s per sample
  Vec2 accel_bias = Vec2(0.08, -0.05);
  double gyro_bias = 0.004;       ///< rad/s (about 0.23 deg/s)
  double accel_bias_walk = 1e-4;  ///< m/s^2 per sqrt(s)
  double gyro_bias_walk = 2e-5;   ///< rad/s per sqrt(s)
  std::uint32_t seed = 11;
};

class ImuSimulator {
 public:
  explicit ImuSimulator(const ImuConfig& cfg = {});

  /// Measurement for the current ego state (CTRV with longitudinal acceleration).
  ImuSample sample(const EgoState& ego, double time_s);

  Vec2 accel_bias() const { return accel_bias_; }
  double gyro_bias() const { return gyro_bias_; }
  const ImuConfig& config() const { return cfg_; }

 private:
  ImuConfig cfg_;
  std::mt19937 rng_;
  Vec2 accel_bias_;
  double gyro_bias_;
};

}  // namespace radar
