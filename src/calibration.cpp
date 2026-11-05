#include "radar/calibration.hpp"

#include <Eigen/Cholesky>
#include <cmath>

namespace radar {

namespace {

Eigen::Matrix3d skew(const Vec3& v) {
  Eigen::Matrix3d s;
  s << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
  return s;
}

/// exp map of se(3) twist (rho, phi), with the exact left Jacobian for translation.
Pose3 se3_exp(const Eigen::Matrix<double, 6, 1>& xi) {
  const Vec3 rho = xi.head<3>(), phi = xi.tail<3>();
  const double th = phi.norm();
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d V = Eigen::Matrix3d::Identity();
  if (th > 1e-12) {
    const Eigen::Matrix3d W = skew(phi);
    R = Eigen::AngleAxisd(th, phi / th).toRotationMatrix();
    V += (1 - std::cos(th)) / (th * th) * W + (th - std::sin(th)) / (th * th * th) * W * W;
  }
  Pose3 T = Pose3::Identity();
  T.linear() = R;
  T.translation() = V * rho;
  return T;
}

double rms(const CameraIntrinsics& K, const std::vector<Correspondence>& data, const Pose3& T) {
  double ss = 0.0;
  for (const auto& c : data) {
    const Vec3 p = T * c.p_radar;
    const Eigen::Vector2d uv(K.fx * p.x() / p.z() + K.cx, K.fy * p.y() / p.z() + K.cy);
    ss += (uv - c.uv).squaredNorm();
  }
  return std::sqrt(ss / static_cast<double>(std::max<std::size_t>(1, data.size())));
}

}  // namespace

CalibrationResult calibrate_extrinsic(const CameraIntrinsics& K,
                                      const std::vector<Correspondence>& data, const Pose3& initial,
                                      int max_iterations) {
  CalibrationResult res;
  Pose3 T = initial;
  res.initial_rms_px = rms(K, data, T);
  double cost = res.initial_rms_px;
  double lambda = 1e-3;

  for (int it = 0; it < max_iterations; ++it) {
    Eigen::Matrix<double, 6, 6> JtJ = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 6, 1> Jtr = Eigen::Matrix<double, 6, 1>::Zero();
    for (const auto& c : data) {
      const Vec3 p = T * c.p_radar;
      if (p.z() <= 0.1) continue;
      const double iz = 1.0 / p.z();
      const Eigen::Vector2d r(K.fx * p.x() * iz + K.cx - c.uv.x(),
                              K.fy * p.y() * iz + K.cy - c.uv.y());
      Eigen::Matrix<double, 2, 3> Jp;
      Jp << K.fx * iz, 0, -K.fx * p.x() * iz * iz, 0, K.fy * iz, -K.fy * p.y() * iz * iz;
      // d(exp(delta) p)/d(delta) = [I, -skew(p)] for delta = (rho, phi).
      Eigen::Matrix<double, 3, 6> Jd;
      Jd << Eigen::Matrix3d::Identity(), -skew(p);
      const Eigen::Matrix<double, 2, 6> J = Jp * Jd;
      JtJ += J.transpose() * J;
      Jtr += J.transpose() * r;
    }
    Eigen::Matrix<double, 6, 6> A = JtJ;
    A.diagonal() += lambda * JtJ.diagonal().cwiseMax(1e-9);
    const Eigen::Matrix<double, 6, 1> delta = -A.ldlt().solve(Jtr);
    const Pose3 candidate = se3_exp(delta) * T;
    const double new_cost = rms(K, data, candidate);
    res.iterations = it + 1;
    if (new_cost < cost) {
      T = candidate;
      lambda = std::max(1e-7, lambda / 10.0);
      const bool small = delta.norm() < 1e-10 || (cost - new_cost) < 1e-9 * cost;
      cost = new_cost;
      if (small) {
        res.converged = true;
        break;
      }
    } else {
      lambda *= 10.0;
      if (lambda > 1e8) {
        res.converged = true;  // no further improvement possible
        break;
      }
    }
  }
  res.T_cam_from_radar = T;
  res.final_rms_px = cost;
  return res;
}

PoseError pose_difference(const Pose3& a, const Pose3& b) {
  const Pose3 d = a.inverse() * b;
  PoseError e;
  e.rotation_rad = Eigen::AngleAxisd(d.linear()).angle();
  e.translation_m = (a.inverse().translation() - b.inverse().translation()).norm();
  return e;
}

}  // namespace radar
