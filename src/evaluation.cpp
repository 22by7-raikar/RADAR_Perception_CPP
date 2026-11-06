#include "radar/evaluation.hpp"

#include <cmath>
#include <limits>

#include "radar/config.hpp"

namespace radar {

EvaluationResult evaluate(const std::vector<Target>& truth,
                          const std::vector<Detection>& detections, const MatchGates& gates) {
  EvaluationResult res;
  res.truth_to_detection.assign(truth.size(), -1);
  std::vector<bool> used(detections.size(), false);

  double sum_r = 0.0, sum_v = 0.0, sum_a = 0.0, sum_e = 0.0;
  for (std::size_t i = 0; i < truth.size(); ++i) {
    const Target& t = truth[i];
    double best = std::numeric_limits<double>::infinity();
    long best_j = -1;
    for (std::size_t j = 0; j < detections.size(); ++j) {
      if (used[j]) continue;
      const double er = (detections[j].range_m - t.range_m) / gates.range_m;
      const double ev = (detections[j].velocity_mps - t.velocity_mps) / gates.velocity_mps;
      if (std::abs(er) > 1.0 || std::abs(ev) > 1.0) continue;
      const double cost = er * er + ev * ev;
      if (cost < best) {
        best = cost;
        best_j = static_cast<long>(j);
      }
    }
    if (best_j < 0) {
      ++res.false_negatives;
      continue;
    }
    used[static_cast<std::size_t>(best_j)] = true;
    res.truth_to_detection[i] = best_j;
    ++res.true_positives;
    const Detection& d = detections[static_cast<std::size_t>(best_j)];
    sum_r += std::abs(d.range_m - t.range_m);
    sum_v += std::abs(d.velocity_mps - t.velocity_mps);
    sum_a += std::abs(d.azimuth_rad - t.azimuth_rad) * 180.0 / kPi;
    sum_e += std::abs(d.elevation_rad - t.elevation_rad) * 180.0 / kPi;
  }
  res.false_positives = detections.size() - res.true_positives;
  if (res.true_positives > 0) {
    const double n = static_cast<double>(res.true_positives);
    res.mean_abs_range_error_m = sum_r / n;
    res.mean_abs_velocity_error_mps = sum_v / n;
    res.mean_abs_azimuth_error_deg = sum_a / n;
    res.mean_abs_elevation_error_deg = sum_e / n;
  }
  return res;
}

}  // namespace radar
