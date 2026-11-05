#pragma once

#include <Eigen/Core>
#include <cmath>

#include "radar/config.hpp"

namespace radar {

using Vec2 = Eigen::Vector2d;
using Mat2 = Eigen::Matrix2d;

/// Wraps an angle to (-pi, pi].
inline double wrap_angle(double a) {
  a = std::fmod(a + kPi, 2.0 * kPi);
  if (a <= 0.0) a += 2.0 * kPi;
  return a - kPi;
}

inline Mat2 rotation(double yaw) {
  const double c = std::cos(yaw), s = std::sin(yaw);
  Mat2 r;
  r << c, -s, s, c;
  return r;
}

/// Planar rigid transform. Maps points from the local frame into the parent
/// frame: p_parent = R(yaw) p_local + t.
struct Pose2D {
  double x = 0.0;
  double y = 0.0;
  double yaw = 0.0;

  Vec2 translation() const { return {x, y}; }
  Mat2 rotation_matrix() const { return rotation(yaw); }
  Vec2 to_parent(const Vec2& local) const { return rotation_matrix() * local + translation(); }
  Vec2 to_local(const Vec2& parent) const {
    return rotation_matrix().transpose() * (parent - translation());
  }
  /// this ∘ child: pose of `child` (given in this frame) expressed in the parent frame.
  Pose2D compose(const Pose2D& child) const {
    const Vec2 t = to_parent(child.translation());
    return {t.x(), t.y(), wrap_angle(yaw + child.yaw)};
  }
};

/// 2D cross product of angular rate w (about +z) with r: w x r.
inline Vec2 cross_z(double w, const Vec2& r) { return {-w * r.y(), w * r.x()}; }

}  // namespace radar
