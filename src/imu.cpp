#include "radar/imu.hpp"

#include <cmath>

namespace radar {

ImuSimulator::ImuSimulator(const ImuConfig& cfg)
    : cfg_(cfg), rng_(cfg.seed), accel_bias_(cfg.accel_bias), gyro_bias_(cfg.gyro_bias) {}

ImuSample ImuSimulator::sample(const EgoState& ego, double time_s) {
  std::normal_distribution<double> n(0.0, 1.0);
  const double dt = 1.0 / cfg_.rate_hz;
  accel_bias_ += cfg_.accel_bias_walk * std::sqrt(dt) * Vec2(n(rng_), n(rng_));
  gyro_bias_ += cfg_.gyro_bias_walk * std::sqrt(dt) * n(rng_);

  // Body-frame acceleration of a vehicle following a CTRV path with
  // longitudinal acceleration: tangential a, centripetal v * w.
  const Vec2 true_accel(ego.accel, ego.speed * ego.yaw_rate);
  ImuSample s;
  s.time_s = time_s;
  s.accel_body = true_accel + accel_bias_ + cfg_.accel_noise_std * Vec2(n(rng_), n(rng_));
  s.gyro_z = ego.yaw_rate + gyro_bias_ + cfg_.gyro_noise_std * n(rng_);
  return s;
}

}  // namespace radar
