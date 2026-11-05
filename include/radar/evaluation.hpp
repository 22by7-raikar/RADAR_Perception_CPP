#pragma once

#include <cstddef>
#include <vector>

#include "radar/simulator.hpp"
#include "radar/types.hpp"

namespace radar {

struct MatchGates {
  double range_m = 2.0;
  double velocity_mps = 2.0;
};

struct EvaluationResult {
  std::size_t true_positives = 0;
  std::size_t false_negatives = 0;
  std::size_t false_positives = 0;
  double mean_abs_range_error_m = 0.0;
  double mean_abs_velocity_error_mps = 0.0;
  double mean_abs_azimuth_error_deg = 0.0;
  /// Index into the detections for each truth target, or -1 if missed.
  std::vector<long> truth_to_detection;
};

/// Greedy nearest-neighbour association of detections to ground truth in
/// (range, velocity), gated per axis.
EvaluationResult evaluate(const std::vector<Target>& truth,
                          const std::vector<Detection>& detections, const MatchGates& gates = {});

}  // namespace radar
