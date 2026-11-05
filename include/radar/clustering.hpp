#pragma once

#include <cstddef>
#include <vector>

#include "radar/geometry.hpp"
#include "radar/types.hpp"

namespace radar {

struct ClusterConfig {
  double eps_position_m = 3.0;    ///< neighbourhood radius in x/y
  double eps_velocity_mps = 2.0;  ///< neighbourhood radius in radial velocity
  std::size_t min_points = 1;     ///< core-point threshold (1 = connected components)
};

/// Object hypothesis built from a group of detections (sensor frame).
struct Cluster {
  std::vector<std::size_t> members;  ///< indices into the detection list
  Vec2 centroid = Vec2::Zero();      ///< RCS-weighted mean position
  double radial_velocity = 0.0;      ///< RCS-weighted mean
  double velocity_spread = 0.0;      ///< max - min radial velocity
  double extent_major = 0.0;         ///< spread along the principal axis (m)
  double extent_minor = 0.0;
  double rcs_dbsm = 0.0;  ///< power sum of member RCS estimates
  double max_snr_db = 0.0;
};

/// DBSCAN over (x, y, v_r) with per-axis neighbourhood radii, i.e. points are
/// neighbours when (dx^2 + dy^2) / eps_p^2 + dv^2 / eps_v^2 <= 1.
///
/// Velocity is part of the metric so that adjacent objects moving differently
/// (a pedestrian next to a parked car) are separated.
class Dbscan {
 public:
  explicit Dbscan(const ClusterConfig& cfg = {}) : cfg_(cfg) {}

  void cluster(const std::vector<Detection>& detections, std::vector<Cluster>& out);

  /// Label per detection from the last call (-1 = noise).
  const std::vector<int>& labels() const { return labels_; }

 private:
  bool neighbours(const Detection& a, const Detection& b) const;

  ClusterConfig cfg_;
  std::vector<int> labels_;
  std::vector<std::size_t> queue_;
  std::vector<std::size_t> scratch_;
};

}  // namespace radar
