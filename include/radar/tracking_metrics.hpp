#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <map>
#include <vector>

#include "radar/scene.hpp"
#include "radar/tracker.hpp"

namespace radar {

struct TrackingEvalConfig {
  double match_distance_m = 1.5;  ///< plus half the object length (centroid vs. centre)
  /// Ground-truth objects slower than this belong to the static environment
  /// and are not scored. Unmatched tracks classified static are not false
  /// positives; matched ones still count in the confusion matrix.
  double dynamic_speed_mps = 0.3;
};

/// CLEAR-MOT style metrics for dynamic objects plus class accuracy.
struct TrackingMetrics {
  std::size_t ground_truth = 0;
  std::size_t matches = 0;
  std::size_t misses = 0;
  std::size_t false_positives = 0;
  std::size_t id_switches = 0;
  double position_error_sum = 0.0;
  double velocity_error_sum = 0.0;
  /// confusion(true class, predicted class), over matched pairs.
  Eigen::Matrix<std::size_t, kNumObjectClasses, kNumObjectClasses> confusion =
      Eigen::Matrix<std::size_t, kNumObjectClasses, kNumObjectClasses>::Zero();

  double mota() const;
  double motp_m() const;
  double mean_velocity_error_mps() const;
  double recall() const;
  /// Fraction of matched pairs whose predicted class equals the true class.
  double class_accuracy() const;
};

class TrackingEvaluator {
 public:
  explicit TrackingEvaluator(const TrackingEvalConfig& cfg = {}) : cfg_(cfg) {}

  void add_frame(const std::vector<GroundTruthObject>& truth, const std::vector<Track>& tracks);
  const TrackingMetrics& metrics() const { return m_; }

 private:
  TrackingEvalConfig cfg_;
  TrackingMetrics m_;
  std::map<int, int> last_track_for_truth_;
};

}  // namespace radar
