#pragma once

#include <Eigen/Core>

#include "radar/geometry.hpp"
#include "radar/imu.hpp"
#include "radar/tracker.hpp"

namespace radar {

struct EgoEstimatorConfig {
  // Continuous-time IMU noise densities and bias random walks.
  double accel_noise = 0.05;  ///< m/s^2 / sqrt(Hz) equivalent
  double gyro_noise = 0.002;  ///< rad/s / sqrt(Hz) equivalent
  double accel_bias_walk = 1e-3;
  double gyro_bias_walk = 1e-4;
  // Initial uncertainty.
  double init_position_std = 0.05;
  double init_velocity_std = 0.3;
  double init_yaw_std = 0.005;
  double init_accel_bias_std = 0.2;
  double init_gyro_bias_std = 0.01;
  /// Radar velocity updates whose innovation exceeds this (2-DoF) are rejected,
  /// e.g. when RANSAC locks onto a group of moving objects.
  double radar_gate_chi2 = 13.8;
  Pose2D sensor_mount;  ///< radar pose in the body frame (extrinsic calibration)
  /// Treat pose as a consider state in radar updates (no pose jumps). Tracking
  /// in the odometry frame needs frame-to-frame consistency more than the last
  /// fraction of a degree of absolute yaw.
  bool consider_pose = true;
};

/// Planar error-state Kalman filter fusing IMU mechanisation with radar
/// Doppler ego-velocity.
///
/// Nominal state: position (world), velocity (body frame), yaw, accelerometer
/// bias (body x/y) and gyro bias. Error state (8): [dp(2), dv_b(2), dyaw,
/// dba(2), dbg]. Keeping velocity in the body frame decouples it from yaw:
/// the radar measures body velocity, so radar updates do not rotate the pose
/// estimate frame by frame.
///
/// The radar measures the sensor's own velocity in the sensor frame,
///   z = R_m^T (R(yaw)^T v + w x r_m),   w = gyro - b_g,
/// so with a mount offset r_m the lateral component also observes yaw rate and
/// therefore the gyro bias.
class EgoEstimator {
 public:
  using ErrorCov = Eigen::Matrix<double, 8, 8>;

  explicit EgoEstimator(const EgoEstimatorConfig& cfg = {});

  void initialize(const Pose2D& pose, const Vec2& velocity_world, double time_s);
  /// Mechanise to the sample's timestamp using the previous sample's rates.
  void propagate(const ImuSample& imu);
  /// Fuses one radar ego-velocity estimate. Returns false if gated out.
  bool update_radar_velocity(const Vec2& sensor_velocity, const Mat2& covariance);

  Pose2D pose() const { return {p_.x(), p_.y(), yaw_}; }
  Vec2 velocity_world() const { return rotation(yaw_) * vb_; }
  Vec2 velocity_body() const { return vb_; }
  double yaw_rate() const { return last_gyro_ - bg_; }
  Vec2 accel_bias() const { return ba_; }
  double gyro_bias() const { return bg_; }
  const ErrorCov& covariance() const { return P_; }
  double time() const { return t_; }
  /// Radar pose and velocity in the world frame, for the object tracker.
  SensorState sensor_state() const;
  bool initialized() const { return initialized_; }

 private:
  EgoEstimatorConfig cfg_;
  Vec2 p_ = Vec2::Zero();
  Vec2 vb_ = Vec2::Zero();
  double yaw_ = 0.0;
  Vec2 ba_ = Vec2::Zero();
  double bg_ = 0.0;
  ErrorCov P_ = ErrorCov::Identity();
  double t_ = 0.0;
  Vec2 last_accel_ = Vec2::Zero();
  double last_gyro_ = 0.0;
  bool have_imu_ = false;
  bool initialized_ = false;
};

/// Radar-only odometry: integrates Doppler ego-velocity, taking yaw rate from
/// the lateral velocity of an offset-mounted sensor (assumes no side slip).
class RadarOdometry {
 public:
  explicit RadarOdometry(const Pose2D& sensor_mount) : mount_(sensor_mount) {}
  void initialize(const Pose2D& pose, double time_s);
  void update(const Vec2& sensor_velocity, double time_s);
  Pose2D pose() const { return pose_; }

 private:
  Pose2D mount_;
  Pose2D pose_;
  double t_ = 0.0;
  double speed_ = 0.0;
  double yaw_rate_ = 0.0;
  bool initialized_ = false;
};

}  // namespace radar
