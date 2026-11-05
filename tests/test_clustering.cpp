#include <gtest/gtest.h>

#include "radar/clustering.hpp"

namespace radar {

namespace {

Detection det(float x, float y, float v, float rcs = 0.0f) {
  Detection d;
  d.x_m = x;
  d.y_m = y;
  d.velocity_mps = v;
  d.rcs_dbsm = rcs;
  return d;
}

}  // namespace

TEST(Dbscan, SeparatesSpatiallyDistinctGroups) {
  std::vector<Detection> d{det(10, 0, 0), det(11, 0.5, 0), det(12, 1, 0),  // object A
                           det(30, 5, 0), det(31, 5, 0)};                  // object B
  Dbscan db;
  std::vector<Cluster> c;
  db.cluster(d, c);
  ASSERT_EQ(c.size(), 2u);
  EXPECT_EQ(c[0].members.size(), 3u);
  EXPECT_EQ(c[1].members.size(), 2u);
  EXPECT_NEAR(c[0].centroid.x(), 11.0, 1e-6);
  EXPECT_NEAR(c[0].extent_major, std::sqrt(4.0 + 1.0), 1e-6);
}

TEST(Dbscan, SeparatesByVelocityAtSamePosition) {
  // A pedestrian standing next to a moving car occupies the same space but
  // has a different radial velocity.
  std::vector<Detection> d{det(20, 0, 8.0), det(21, 0.5, 8.1), det(20.5, 0.8, 0.2)};
  Dbscan db;
  std::vector<Cluster> c;
  db.cluster(d, c);
  ASSERT_EQ(c.size(), 2u);
  EXPECT_EQ(db.labels()[0], db.labels()[1]);
  EXPECT_NE(db.labels()[0], db.labels()[2]);
}

TEST(Dbscan, ChainsPointsAlongExtendedObject) {
  std::vector<Detection> d;
  for (int i = 0; i < 5; ++i) d.push_back(det(50.0f + 2.0f * static_cast<float>(i), 3.0f, -5.0f));
  Dbscan db;
  std::vector<Cluster> c;
  db.cluster(d, c);
  ASSERT_EQ(c.size(), 1u);
  EXPECT_NEAR(c[0].extent_major, 8.0, 1e-6);
  EXPECT_NEAR(c[0].extent_minor, 0.0, 1e-6);
}

TEST(Dbscan, RcsWeightedCentroidAndPowerSum) {
  std::vector<Detection> d{det(0, 0, 0, 10.0f), det(1, 0, 0, 0.0f)};
  Dbscan db;
  std::vector<Cluster> c;
  db.cluster(d, c);
  ASSERT_EQ(c.size(), 1u);
  EXPECT_NEAR(c[0].centroid.x(), 1.0 / 11.0, 1e-6);  // weights 10 and 1
  EXPECT_NEAR(c[0].rcs_dbsm, 10.0 * std::log10(11.0), 1e-6);
}

TEST(Dbscan, MinPointsMarksNoise) {
  ClusterConfig cfg;
  cfg.min_points = 3;
  std::vector<Detection> d{det(0, 0, 0), det(0.5, 0, 0), det(1, 0, 0), det(40, 40, 0)};
  Dbscan db(cfg);
  std::vector<Cluster> c;
  db.cluster(d, c);
  ASSERT_EQ(c.size(), 1u);
  EXPECT_EQ(db.labels()[3], -1);
}

TEST(Dbscan, EmptyInput) {
  Dbscan db;
  std::vector<Cluster> c{Cluster{}};
  db.cluster({}, c);
  EXPECT_TRUE(c.empty());
}

}  // namespace radar
