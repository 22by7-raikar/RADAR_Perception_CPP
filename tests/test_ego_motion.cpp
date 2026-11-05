#include <gtest/gtest.h>

#include <Eigen/LU>
#include <cmath>
#include <random>

#include "radar/ego_estimator.hpp"
#include "radar/ego_velocity.hpp"
#include "radar/imu.hpp"

namespace radar {

namespace {

constexpr double kDeg = kPi / 180.0;

Detection static_point(double az, const Vec2& v_sensor, double noise, std::mt19937& rng) {
  std::normal_distribution<double> n(0.0, noise);
  Detection d;
  d.azimuth_rad = static_cast<float>(az);
  d.velocity_mps =
      static_cast<float>(-(std::cos(az) * v_sensor.x() + std::sin(az) * v_sensor.y()) + n(rng));
  return d;
}

}  // namespace

TEST(EgoVelocity, RecoversSensorVelocityDespiteMovingObjects) {
  std::mt19937 rng(1);
  const Vec2 v(12.0, -0.4);
  std::vector<Detection> d;
  for (int i = 0; i < 30; ++i) d.push_back(static_point((-60.0 + 4.0 * i) * kDeg, v, 0.05, rng));
  // A third of the points are moving objects.
  for (int i = 0; i < 12; ++i) {
    Detection m = static_point((-50.0 + 8.0 * i) * kDeg, v, 0.05, rng);
    m.velocity_mps += static_cast<float>(3.0 + i);
    d.push_back(m);
  }
  EgoVelocityEstimator est;
  const auto e = est.estimate(d);
  ASSERT_TRUE(e.valid);
  EXPECT_NEAR(e.velocity.x(), v.x(), 0.03);
  EXPECT_NEAR(e.velocity.y(), v.y(), 0.05);
  EXPECT_EQ(e.inliers, 30u);
  for (int i = 30; i < 42; ++i) EXPECT_EQ(e.inlier_mask[static_cast<std::size_t>(i)], 0);
  EXPECT_GT(e.covariance.determinant(), 0.0);
}

TEST(EgoVelocity, InvalidWithTooFewOrDegeneratePoints) {
  std::mt19937 rng(2);
  EgoVelocityEstimator est;
  std::vector<Detection> few;
  for (int i = 0; i < 3; ++i) few.push_back(static_point(0.1 * i, {10, 0}, 0.0, rng));
  EXPECT_FALSE(est.estimate(few).valid);

  std::vector<Detection> same_angle;
  for (int i = 0; i < 20; ++i) same_angle.push_back(static_point(0.2, {10, 0}, 0.0, rng));
  EXPECT_FALSE(est.estimate(same_angle).valid);
}

TEST(Imu, NoiseFreeSamplesMatchKinematics) {
  ImuConfig cfg;
  cfg.accel_noise_std = cfg.gyro_noise_std = 0.0;
  cfg.accel_bias = Vec2::Zero();
  cfg.gyro_bias = 0.0;
  cfg.accel_bias_walk = cfg.gyro_bias_walk = 0.0;
  ImuSimulator imu(cfg);
  EgoState e;
  e.speed = 10.0;
  e.yaw_rate = 0.2;
  e.accel = 1.5;
  const ImuSample s = imu.sample(e, 0.0);
  EXPECT_NEAR(s.accel_body.x(), 1.5, 1e-12);
  EXPECT_NEAR(s.accel_body.y(), 2.0, 1e-12);  // centripetal v * w
  EXPECT_NEAR(s.gyro_z, 0.2, 1e-12);
}

TEST(EgoEstimator, NoiseFreeMechanisationFollowsCircle) {
  // Constant speed and turn rate: the filter must reproduce the CTRV arc.
  EgoEstimatorConfig cfg;
  EgoEstimator f(cfg);
  EgoState e;
  e.speed = 10.0;
  e.yaw_rate = 0.1;
  f.initialize(e.pose, e.velocity(), 0.0);
  Vec2 p = Vec2::Zero();
  double yaw = 0.0;
  for (int k = 0; k <= 1000; ++k) {
    const double t = 0.01 * k;
    f.propagate({t, Vec2(0.0, e.speed * e.yaw_rate), e.yaw_rate});
    if (k < 1000) propagate_ctrv(p, yaw, e.speed, e.yaw_rate, 0.01);
  }
  EXPECT_LT((f.pose().translation() - p).norm(), 0.05);
  EXPECT_NEAR(f.pose().yaw, yaw, 1e-6);
}

TEST(EgoEstimator, RadarLeverArmMakesGyroBiasObservable) {
  // Note: at constant speed in a straight line, gyro bias, lateral velocity and
  // lateral accelerometer bias are not separable (b_g * v_x and b_a,y cancel).
  // A changing speed makes them observable, as in normal driving.
  EgoEstimatorConfig cfg;
  cfg.sensor_mount = {3.6, 0.0, 0.0};
  EgoEstimator f(cfg);
  const double accel = 0.5, gyro_bias = 0.01;
  double speed = 8.0;
  f.initialize({}, {speed, 0.0}, 0.0);
  std::mt19937 rng(4);
  std::normal_distribution<double> n(0.0, 1.0);
  for (int k = 0; k <= 3000; ++k) {
    const double t = 0.01 * k;
    f.propagate({t, Vec2(accel, 0.0), gyro_bias + 0.002 * n(rng)});
    if (k % 5 == 0) {
      // Straight driving: the radar sees (v, 0); any lateral motion would be w * r_x.
      const Vec2 z(speed + 0.05 * n(rng), 0.05 * n(rng));
      f.update_radar_velocity(z, 0.05 * 0.05 * Mat2::Identity());
    }
    speed += accel * 0.01;
  }
  EXPECT_NEAR(f.gyro_bias(), gyro_bias, 0.002);
  EXPECT_NEAR(f.velocity_body().x(), speed, 0.1);
  EXPECT_NEAR(f.velocity_body().y(), 0.0, 0.1);
}

TEST(EgoEstimator, ConsiderUpdateDoesNotMovePose) {
  EgoEstimatorConfig cfg;
  EgoEstimator f(cfg);
  f.initialize({5.0, 2.0, 0.3}, 10.0 * Vec2(std::cos(0.3), std::sin(0.3)), 0.0);
  for (int k = 0; k <= 50; ++k) f.propagate({0.01 * k, Vec2::Zero(), 0.0});
  const Pose2D before = f.pose();
  EXPECT_TRUE(f.update_radar_velocity({9.8, 0.1}, 0.1 * Mat2::Identity()));
  EXPECT_DOUBLE_EQ(f.pose().x, before.x);
  EXPECT_DOUBLE_EQ(f.pose().y, before.y);
  EXPECT_DOUBLE_EQ(f.pose().yaw, before.yaw);
}

TEST(EgoEstimator, GatesOutlierVelocity) {
  EgoEstimatorConfig cfg;
  EgoEstimator f(cfg);
  f.initialize({}, {10.0, 0.0}, 0.0);
  f.propagate({0.0, Vec2::Zero(), 0.0});
  f.propagate({0.01, Vec2::Zero(), 0.0});
  EXPECT_FALSE(f.update_radar_velocity({25.0, 0.0}, 0.01 * Mat2::Identity()));
}

TEST(RadarOdometry, YawRateFromOffsetMount) {
  const Pose2D mount{3.6, 0.0, 0.0};
  RadarOdometry odo(mount);
  odo.initialize({}, 0.0);
  const double v = 10.0, w = 0.2;
  Vec2 p = Vec2::Zero();
  double yaw = 0.0;
  for (int k = 1; k <= 100; ++k) {
    // Sensor velocity in the sensor frame for a CTRV vehicle.
    odo.update({v, w * mount.x}, 0.05 * k);
    propagate_ctrv(p, yaw, v, w, 0.05);
  }
  // First interval integrates zero velocity, so allow one step of lag.
  EXPECT_LT((odo.pose().translation() - p).norm(), v * 0.05 + 0.1);
  EXPECT_NEAR(odo.pose().yaw, yaw, w * 0.05 + 1e-6);
}

TEST(Geometry, PoseInverse) {
  const Pose2D p{4.0, -1.0, 0.7};
  const Pose2D id = p.compose(p.inverse());
  EXPECT_NEAR(id.x, 0.0, 1e-12);
  EXPECT_NEAR(id.y, 0.0, 1e-12);
  EXPECT_NEAR(id.yaw, 0.0, 1e-12);
}

}  // namespace radar
