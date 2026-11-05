#include "radar/tracking_metrics.hpp"

#include <algorithm>

#include "radar/assignment.hpp"

namespace radar {

double TrackingMetrics::mota() const {
  if (ground_truth == 0) return 1.0;
  return 1.0 - static_cast<double>(misses + false_positives + id_switches) /
                   static_cast<double>(ground_truth);
}

double TrackingMetrics::motp_m() const {
  return matches ? position_error_sum / static_cast<double>(matches) : 0.0;
}

double TrackingMetrics::mean_velocity_error_mps() const {
  return matches ? velocity_error_sum / static_cast<double>(matches) : 0.0;
}

double TrackingMetrics::recall() const {
  return ground_truth ? static_cast<double>(matches) / static_cast<double>(ground_truth) : 0.0;
}

double TrackingMetrics::class_accuracy() const {
  const std::size_t total = confusion.sum();
  return total ? static_cast<double>(confusion.trace()) / static_cast<double>(total) : 0.0;
}

void TrackingEvaluator::add_frame(const std::vector<GroundTruthObject>& truth,
                                  const std::vector<Track>& tracks, const Pose2D& track_to_truth) {
  std::vector<const GroundTruthObject*> gt;
  for (const auto& g : truth) {
    if (g.velocity_world.norm() >= cfg_.dynamic_speed_mps) gt.push_back(&g);
  }
  // Every confirmed track may match a dynamic object (a misclassified track
  // is a class error, not a miss); unmatched tracks are false positives only
  // when they claim to be dynamic.
  std::vector<const Track*> trk;
  std::vector<Vec2> trk_pos, trk_vel;
  const Mat2 rot = track_to_truth.rotation_matrix();
  for (const auto& t : tracks) {
    if (!t.confirmed()) continue;
    trk.push_back(&t);
    trk_pos.push_back(track_to_truth.to_parent(t.position()));
    trk_vel.push_back(rot * t.velocity());
  }

  constexpr double kForbidden = 1e9;
  Eigen::MatrixXd cost = Eigen::MatrixXd::Constant(
      static_cast<Eigen::Index>(gt.size()), static_cast<Eigen::Index>(trk.size()), kForbidden);
  for (std::size_t i = 0; i < gt.size(); ++i) {
    for (std::size_t j = 0; j < trk.size(); ++j) {
      const double d = (gt[i]->position_world - trk_pos[j]).norm();
      // Gate grows with object size: radar sees the near surface, not the centre.
      const double gate = cfg_.match_distance_m + 0.5 * std::max(gt[i]->length, gt[i]->width);
      if (d <= gate) {
        cost(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = d;
      }
    }
  }
  const auto match = solve_assignment(cost, kForbidden);

  m_.ground_truth += gt.size();
  for (std::size_t i = 0; i < gt.size(); ++i) {
    if (match[i] < 0) {
      ++m_.misses;
      continue;
    }
    const auto j = static_cast<std::size_t>(match[i]);
    const Track& t = *trk[j];
    ++m_.matches;
    m_.position_error_sum += (gt[i]->position_world - trk_pos[j]).norm();
    m_.velocity_error_sum += (gt[i]->velocity_world - trk_vel[j]).norm();
    m_.confusion(static_cast<Eigen::Index>(gt[i]->cls), static_cast<Eigen::Index>(t.cls)) += 1;
    auto it = last_track_for_truth_.find(gt[i]->id);
    if (it != last_track_for_truth_.end() && it->second != t.id) ++m_.id_switches;
    last_track_for_truth_[gt[i]->id] = t.id;
  }
  std::vector<bool> taken(trk.size(), false);
  for (int j : match) {
    if (j >= 0) taken[static_cast<std::size_t>(j)] = true;
  }
  for (std::size_t j = 0; j < trk.size(); ++j) {
    if (!taken[j] && trk[j]->cls != ObjectClass::kStatic) ++m_.false_positives;
  }
}

}  // namespace radar
