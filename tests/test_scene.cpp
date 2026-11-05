#include <gtest/gtest.h>

#include <cmath>

#include "radar/scene.hpp"

namespace radar {

constexpr double kDeg = kPi / 180.0;

TEST(Geometry, PoseRoundTripAndCompose) {
  const Pose2D p{3.0, -2.0, 30.0 * kDeg};
  const Vec2 local(5.0, 1.0);
  EXPECT_TRUE(p.to_local(p.to_parent(local)).isApprox(local, 1e-12));
  const Pose2D child{1.0, 0.0, 15.0 * kDeg};
  const Pose2D c = p.compose(child);
  EXPECT_TRUE(c.to_parent(local).isApprox(p.to_parent(child.to_parent(local)), 1e-12));
  EXPECT_NEAR(wrap_angle(3.5 * kPi), -0.5 * kPi, 1e-12);
}

TEST(Scene, CtrvStraightAndTurning) {
  Vec2 p(0, 0);
  double h = 0.0;
  propagate_ctrv(p, h, 10.0, 0.0, 2.0);
  EXPECT_TRUE(p.isApprox(Vec2(20.0, 0.0)));

  // Quarter circle of radius 10.
  p = Vec2(0, 0);
  h = 0.0;
  propagate_ctrv(p, h, 10.0 * kPi / 2.0, 1.0, kPi / 2.0 / 1.0);
  EXPECT_NEAR(h, kPi / 2.0, 1e-9);
}

TEST(Scene, StaticReflectorSeenFromMovingEgo) {
  // Radial velocity of a static point is -u . v_ego.
  EgoState ego;
  ego.speed = 10.0;
  Scene s({}, ego);
  s.add(SceneObject::pole(1, {50.0, 50.0}));
  const SceneFrame f = s.snapshot();
  ASSERT_EQ(f.targets.size(), 1u);
  EXPECT_NEAR(f.targets[0].range_m, std::hypot(50.0, 50.0), 1e-9);
  EXPECT_NEAR(f.targets[0].azimuth_rad, 45.0 * kDeg, 1e-12);
  EXPECT_NEAR(f.targets[0].velocity_mps, -10.0 * std::cos(45.0 * kDeg), 1e-9);
}

TEST(Scene, RotationOfEgoDoesNotChangeRadialVelocity) {
  EgoState ego;
  ego.yaw_rate = 0.5;  // turning on the spot
  Scene s({}, ego);
  s.add(SceneObject::pole(1, {30.0, 10.0}));
  EXPECT_NEAR(s.snapshot().targets[0].velocity_mps, 0.0, 1e-12);
}

TEST(Scene, MountOffsetAddsLeverArmVelocity) {
  SceneConfig cfg;
  cfg.sensor_mount = {2.0, 0.0, 0.0};  // 2 m ahead of the rotation centre
  EgoState ego;
  ego.yaw_rate = 0.5;
  Scene s(cfg, ego);
  const SceneFrame f = s.snapshot();
  EXPECT_TRUE(f.sensor_velocity.isApprox(Vec2(0.0, 1.0)));  // w x r
}

TEST(Scene, VehicleShowsOnlyFacingScatterers) {
  Scene s;
  s.add(SceneObject::vehicle(1, {40.0, 0.0}, 0.0, 0.0));  // rear faces the sensor
  const SceneFrame f = s.snapshot();
  ASSERT_EQ(f.objects.size(), 1u);
  // Rear centre and the two rear corners are visible; the front is not.
  EXPECT_GE(f.objects[0].num_scatterers, 3u);
  EXPECT_LT(f.objects[0].num_scatterers, 10u);
  for (const auto& t : f.targets) EXPECT_LT(t.range_m, 40.0);
}

TEST(Scene, FieldOfViewAndRangeLimits) {
  Scene s;
  s.add(SceneObject::pole(1, {0.0, 30.0}));   // 90 deg off boresight
  s.add(SceneObject::pole(2, {300.0, 0.0}));  // beyond max range
  s.add(SceneObject::pole(3, {30.0, 0.0}));
  const SceneFrame f = s.snapshot();
  ASSERT_EQ(f.objects.size(), 1u);
  EXPECT_EQ(f.objects[0].id, 3);
}

TEST(Scene, PedestrianLimbsSpreadDoppler) {
  Scene s;
  s.add(SceneObject::pedestrian(1, {20.0, 0.0}, 0.0, 1.5));  // walking away
  s.step(0.13);
  const SceneFrame f = s.snapshot();
  ASSERT_EQ(f.targets.size(), 3u);
  double lo = 1e9, hi = -1e9;
  for (const auto& t : f.targets) {
    lo = std::min(lo, t.velocity_mps);
    hi = std::max(hi, t.velocity_mps);
  }
  EXPECT_GT(hi - lo, 1.0);
}

}  // namespace radar
