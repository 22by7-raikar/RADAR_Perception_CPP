#include "radar/ego_estimator.hpp"

#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <cmath>

namespace radar {

namespace {

/// d/dyaw of R(yaw).
Mat2 rotation_derivative(double yaw) {
  const double c = std::cos(yaw), s = std::sin(yaw);
  Mat2 d;
  d << -s, -c, c, -s;
  return d;
}

}  // namespace

EgoEstimator::EgoEstimator(const EgoEstimatorConfig& cfg) : cfg_(cfg) {}

void EgoEstimator::initialize(const Pose2D& pose, const Vec2& velocity_world, double time_s) {
  p_ = pose.translation();
  yaw_ = pose.yaw;
  vb_ = rotation(yaw_).transpose() * velocity_world;
  ba_.setZero();
  bg_ = 0.0;
  t_ = time_s;
  P_.setZero();
  P_.diagonal() << cfg_.init_position_std, cfg_.init_position_std, cfg_.init_velocity_std,
      cfg_.init_velocity_std, cfg_.init_yaw_std, cfg_.init_accel_bias_std, cfg_.init_accel_bias_std,
      cfg_.init_gyro_bias_std;
  P_ = P_.cwiseProduct(P_).eval();
  have_imu_ = false;
  initialized_ = true;
}

void EgoEstimator::propagate(const ImuSample& imu) {
  if (!have_imu_) {
    last_accel_ = imu.accel_body;
    last_gyro_ = imu.gyro_z;
    t_ = imu.time_s;
    have_imu_ = true;
    return;
  }
  const double dt = imu.time_s - t_;
  if (dt > 0.0) {
    const double w = last_gyro_ - bg_;
    const Vec2 a = last_accel_ - ba_;
    const double yaw_mid = yaw_ + 0.5 * w * dt;
    const Mat2 R = rotation(yaw_mid);

    // Body-frame kinematics: dv_b/dt = a - w x v_b, dp/dt = R v_b, dyaw/dt = w.
    ErrorCov F = ErrorCov::Identity();
    F.block<2, 2>(0, 2) = R * dt;
    F.block<2, 1>(0, 4) = rotation_derivative(yaw_mid) * vb_ * dt;
    Mat2 Wx;
    Wx << 0.0, w, -w, 0.0;  // d(-w x v)/dv
    F.block<2, 2>(2, 2) += Wx * dt;
    F.block<2, 2>(2, 5) = -Mat2::Identity() * dt;
    F.block<2, 1>(2, 7) = Vec2(-vb_.y(), vb_.x()) * dt;  // d(-w x v)/db_g with w = gyro - b_g
    F(4, 7) = -dt;
    ErrorCov Q = ErrorCov::Zero();
    const double qa = cfg_.accel_noise * dt, qg = cfg_.gyro_noise * dt;
    Q(2, 2) = Q(3, 3) = qa * qa;
    Q(4, 4) = qg * qg;
    Q(5, 5) = Q(6, 6) = cfg_.accel_bias_walk * cfg_.accel_bias_walk * dt;
    Q(7, 7) = cfg_.gyro_bias_walk * cfg_.gyro_bias_walk * dt;

    const Vec2 vb_next = vb_ + (a - cross_z(w, vb_)) * dt;
    p_ += R * (0.5 * (vb_ + vb_next)) * dt;
    vb_ = vb_next;
    yaw_ = wrap_angle(yaw_ + w * dt);
    P_ = F * P_ * F.transpose() + Q;
    t_ = imu.time_s;
  }
  last_accel_ = imu.accel_body;
  last_gyro_ = imu.gyro_z;
}

bool EgoEstimator::update_radar_velocity(const Vec2& z, const Mat2& Rz) {
  const Mat2 Rm_t = cfg_.sensor_mount.rotation_matrix().transpose();
  const Vec2 r = cfg_.sensor_mount.translation();
  const double w = last_gyro_ - bg_;
  const Vec2 z_pred = Rm_t * (vb_ + cross_z(w, r));

  // The radar observes body velocity and (through the lever arm) yaw rate;
  // pose is left to the mechanisation, so updates never make it jump.
  Eigen::Matrix<double, 2, 8> H = Eigen::Matrix<double, 2, 8>::Zero();
  H.block<2, 2>(0, 2) = Rm_t;
  H.block<2, 1>(0, 7) = Rm_t * Vec2(r.y(), -r.x());  // d(w x r)/d(b_g)

  const Vec2 y = z - z_pred;
  const Mat2 S = H * P_ * H.transpose() + Rz;
  const Eigen::LDLT<Mat2> ldlt(S);
  if (y.dot(ldlt.solve(y)) > cfg_.radar_gate_chi2) return false;

  Eigen::Matrix<double, 8, 2> K = P_ * H.transpose() * S.inverse();
  if (cfg_.consider_pose) {
    // Schmidt (consider) update: pose is not corrected by velocity updates, so
    // it evolves smoothly from the corrected velocity and gyro bias. The Joseph
    // form below keeps P exact for this suboptimal gain.
    K.row(0).setZero();
    K.row(1).setZero();
    K.row(4).setZero();
  }
  const Eigen::Matrix<double, 8, 1> dx = K * y;
  p_ += dx.segment<2>(0);
  vb_ += dx.segment<2>(2);
  yaw_ = wrap_angle(yaw_ + dx(4));
  ba_ += dx.segment<2>(5);
  bg_ += dx(7);

  const ErrorCov I_KH = ErrorCov::Identity() - K * H;
  P_ = I_KH * P_ * I_KH.transpose() + K * Rz * K.transpose();
  P_ = 0.5 * (P_ + P_.transpose());
  return true;
}

SensorState EgoEstimator::sensor_state() const {
  const Pose2D body = pose();
  SensorState s;
  s.pose = body.compose(cfg_.sensor_mount);
  s.velocity = velocity_world() +
               cross_z(yaw_rate(), body.rotation_matrix() * cfg_.sensor_mount.translation());
  return s;
}

void RadarOdometry::initialize(const Pose2D& pose, double time_s) {
  pose_ = pose;
  t_ = time_s;
  speed_ = 0.0;
  yaw_rate_ = 0.0;
  initialized_ = true;
}

void RadarOdometry::update(const Vec2& sensor_velocity, double time_s) {
  if (!initialized_) return;
  // Integrate with the previous estimate, then take the new one.
  Vec2 p = pose_.translation();
  double yaw = pose_.yaw;
  propagate_ctrv(p, yaw, speed_, yaw_rate_, time_s - t_);
  pose_ = {p.x(), p.y(), yaw};
  t_ = time_s;

  // Sensor velocity in the body frame; with no side slip, v_body = (u, 0), so
  //   v_s = (u - w r_y, w r_x).
  const Vec2 vb = mount_.rotation_matrix() * sensor_velocity;
  const Vec2 r = mount_.translation();
  if (std::abs(r.x()) > 0.1) yaw_rate_ = vb.y() / r.x();
  speed_ = vb.x() + yaw_rate_ * r.y();
}

}  // namespace radar
