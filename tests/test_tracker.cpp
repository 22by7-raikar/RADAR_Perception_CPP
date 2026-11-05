#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "radar/tracker.hpp"

namespace radar {

namespace {

/// Cluster as the radar would report it for a point at `p` (world) moving at
/// `v`, seen from a sensor at the origin with optional ego velocity.
Cluster observe(const Vec2& p, const Vec2& v, const SensorState& s, std::mt19937& rng,
                double pos_noise = 0.3, double vr_noise = 0.1, double rcs = 15.0,
                double spread = 0.0) {
  std::normal_distribution<double> n(0.0, 1.0);
  Cluster c;
  c.centroid = s.pose.to_local(p) + Vec2(pos_noise * n(rng), pos_noise * n(rng));
  const Vec2 d = p - s.pose.translation();
  c.radial_velocity = d.dot(v - s.velocity) / d.norm() + vr_noise * n(rng);
  c.rcs_dbsm = rcs;
  c.extent_major = rcs > 8.0 ? 3.0 : 0.3;
  c.velocity_spread = spread;
  return c;
}

}  // namespace

TEST(Tracker, ConfirmsAfterThreeHitsAndConvergesOnRadialMotion) {
  MultiTargetTracker trk;
  std::mt19937 rng(1);
  SensorState s;
  Vec2 p(60.0, 2.0);
  const Vec2 v(-10.0, 0.0);
  for (int k = 0; k < 40; ++k) {
    trk.update(0.05 * k, {observe(p, v, s, rng)}, s);
    if (k == 1) {
      EXPECT_TRUE(trk.confirmed_tracks().empty());
    }
    if (k == 2) {
      EXPECT_EQ(trk.confirmed_tracks().size(), 1u);
    }
    p += 0.05 * v;
  }
  const auto c = trk.confirmed_tracks();
  ASSERT_EQ(c.size(), 1u);
  EXPECT_LT((c[0].position() - p + 0.05 * v).norm(), 0.5);
  EXPECT_LT((c[0].velocity() - v).norm(), 0.3);
  EXPECT_EQ(c[0].cls, ObjectClass::kVehicle);
}

TEST(Tracker, RecoversTangentialVelocityFromPositionHistory) {
  // Crossing target: radial velocity ~0, motion only visible in position.
  MultiTargetTracker trk;
  std::mt19937 rng(2);
  SensorState s;
  Vec2 p(40.0, -15.0);
  const Vec2 v(0.0, 8.0);
  for (int k = 0; k < 80; ++k) {
    trk.update(0.05 * k, {observe(p, v, s, rng)}, s);
    p += 0.05 * v;
  }
  const auto c = trk.confirmed_tracks();
  ASSERT_EQ(c.size(), 1u);
  EXPECT_LT((c[0].velocity() - v).norm(), 1.0);
}

TEST(Tracker, CompensatesEgoMotion) {
  // Static object seen from a moving sensor must be tracked as static.
  MultiTargetTracker trk;
  std::mt19937 rng(3);
  const Vec2 pole(50.0, 5.0);
  for (int k = 0; k < 60; ++k) {
    SensorState s;
    s.pose = {0.6 * k, 0.0, 0.0};
    s.velocity = Vec2(12.0, 0.0);
    trk.update(0.05 * k, {observe(pole, Vec2::Zero(), s, rng, 0.3, 0.1, 5.0)}, s);
  }
  const auto c = trk.confirmed_tracks();
  ASSERT_EQ(c.size(), 1u);
  EXPECT_LT(c[0].velocity().norm(), 0.5);
  EXPECT_LT((c[0].position() - pole).norm(), 0.5);
  EXPECT_EQ(c[0].cls, ObjectClass::kStatic);
}

TEST(Tracker, DeletesTrackAfterMisses) {
  TrackerConfig cfg;
  MultiTargetTracker trk(cfg);
  std::mt19937 rng(4);
  SensorState s;
  int k = 0;
  for (; k < 5; ++k) trk.update(0.05 * k, {observe({30, 0}, {0, 0}, s, rng)}, s);
  ASSERT_EQ(trk.confirmed_tracks().size(), 1u);
  for (int m = 0; m <= cfg.max_confirmed_misses; ++m, ++k) trk.update(0.05 * k, {}, s);
  EXPECT_TRUE(trk.tracks().empty());
}

TEST(Tracker, KeepsIdentitiesOfTwoCrossingTargets) {
  MultiTargetTracker trk;
  std::mt19937 rng(5);
  SensorState s;
  Vec2 a(30.0, -10.0), b(30.0, 10.0);
  const Vec2 va(0.0, 6.0), vb(0.0, -6.0);
  int id_a = -1, id_b = -1;
  for (int k = 0; k < 70; ++k) {
    // Different radial velocities keep them separable in the Doppler domain.
    trk.update(0.05 * k,
               {observe(a, va + Vec2(4, 0), s, rng), observe(b, vb + Vec2(-4, 0), s, rng)}, s);
    a += 0.05 * (va + Vec2(4, 0));
    b += 0.05 * (vb + Vec2(-4, 0));
    if (k == 10) {
      for (const auto& t : trk.confirmed_tracks()) (t.position().y() < 0 ? id_a : id_b) = t.id;
    }
  }
  ASSERT_NE(id_a, -1);
  ASSERT_NE(id_b, -1);
  for (const auto& t : trk.confirmed_tracks()) {
    if (t.id == id_a) {
      EXPECT_GT(t.velocity().x(), 2.0);
    }
    if (t.id == id_b) {
      EXPECT_LT(t.velocity().x(), -2.0);
    }
  }
}

TEST(Tracker, SecondClusterOnSameVehicleDoesNotSpawnTrack) {
  MultiTargetTracker trk;
  std::mt19937 rng(6);
  SensorState s;
  Vec2 p(50.0, 0.0);
  const Vec2 v(-8.0, 0.0);
  for (int k = 0; k < 30; ++k) {
    // Front and rear of a vehicle sometimes resolve as two clusters 4 m apart.
    std::vector<Cluster> cl{observe(p, v, s, rng)};
    if (k % 3 == 0) cl.push_back(observe(p + Vec2(4.0, 0.0), v, s, rng));
    trk.update(0.05 * k, cl, s);
    p += 0.05 * v;
  }
  EXPECT_EQ(trk.confirmed_tracks().size(), 1u);
}

TEST(Tracker, ClassifiesSmallSpreadingTargetAsPedestrian) {
  MultiTargetTracker trk;
  std::mt19937 rng(7);
  SensorState s;
  Vec2 p(25.0, 3.0);
  const Vec2 v(-1.4, 0.0);
  for (int k = 0; k < 40; ++k) {
    trk.update(0.05 * k, {observe(p, v, s, rng, 0.2, 0.1, -2.0, 1.5)}, s);
    p += 0.05 * v;
  }
  const auto c = trk.confirmed_tracks();
  ASSERT_EQ(c.size(), 1u);
  EXPECT_EQ(c[0].cls, ObjectClass::kPedestrian);
  EXPECT_GT(c[0].class_confidence, 0.5);
}

}  // namespace radar
