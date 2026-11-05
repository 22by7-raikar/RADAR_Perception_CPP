#include "radar/fusion.hpp"

#include <algorithm>
#include <cmath>

#include "radar/assignment.hpp"

namespace radar {

RadarCameraFusion::RadarCameraFusion(const CameraModel& camera, const FusionConfig& cfg)
    : camera_(camera), cfg_(cfg) {}

std::vector<FusedObject> RadarCameraFusion::fuse(const std::vector<Track>& tracks,
                                                 double radar_time_s, const SensorState& sensor,
                                                 double camera_time_s,
                                                 const std::vector<CameraDetection>& dets) {
  std::vector<FusedObject> out;
  const double dt = cfg_.compensate_time ? camera_time_s - radar_time_s : 0.0;

  for (const Track& t : tracks) {
    if (!t.confirmed()) continue;
    FusedObject f;
    f.track_id = t.id;
    f.velocity_world = t.velocity();
    const Vec2 p_world = t.position() + t.velocity() * dt;
    f.position_sensor = sensor.pose.to_local(p_world);

    // Box: class prior for height/width, radar extent for length when larger.
    // Size prior: a "static" radar track that the camera can see is most
    // likely a slow pedestrian, so use a pedestrian-sized box for it.
    const ObjectSize prior =
        nominal_size(t.cls == ObjectClass::kStatic ? ObjectClass::kPedestrian : t.cls);
    const double speed = t.velocity().norm();
    const double heading_world =
        speed > 1.0 ? std::atan2(t.velocity().y(), t.velocity().x()) : sensor.pose.yaw;
    const double heading = wrap_angle(heading_world - sensor.pose.yaw);
    const double length = std::max(prior.length, t.extent_m);
    // The radar sees the near surface; shift the centre back by half a length
    // along the line of sight so the box covers the whole object.
    const Vec2 los = f.position_sensor.normalized();
    const Vec2 centre = f.position_sensor + 0.5 * prior.width * los;
    if (!camera_.project_box(centre, heading, length, prior.width, prior.height,
                             -cfg_.radar_height_m, f.projected)) {
      continue;
    }
    f.cls = t.cls;
    f.class_confidence = t.class_confidence;
    out.push_back(f);
  }

  constexpr double kForbidden = 1e9;
  Eigen::MatrixXd cost = Eigen::MatrixXd::Constant(
      static_cast<Eigen::Index>(out.size()), static_cast<Eigen::Index>(dets.size()), kForbidden);
  for (std::size_t i = 0; i < out.size(); ++i) {
    for (std::size_t j = 0; j < dets.size(); ++j) {
      const double iou = out[i].projected.iou(dets[j].box);
      if (iou >= cfg_.min_iou) {
        cost(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = 1.0 - iou;
      }
    }
  }
  const auto match = solve_assignment(cost, kForbidden);

  for (std::size_t i = 0; i < out.size(); ++i) {
    FusedObject& f = out[i];
    auto& ev = camera_evidence_[f.track_id];
    for (double& e : ev) e *= cfg_.evidence_decay;
    if (match[i] >= 0) {
      const auto j = static_cast<std::size_t>(match[i]);
      f.camera_index = match[i];
      f.iou = f.projected.iou(dets[j].box);
      ev[static_cast<std::size_t>(dets[j].cls)] += dets[j].score;
    }
    // Radar class probability (confidence on its class, the rest spread
    // evenly), plus normalised camera evidence.
    double ev_sum = 1.0;
    for (double e : ev) ev_sum += e;
    std::array<double, kNumObjectClasses> score{};
    const double p_own = f.cls == ObjectClass::kUnknown ? 0.25 : f.class_confidence;
    for (std::size_t k = 0; k < score.size(); ++k) {
      const double p_radar =
          (static_cast<ObjectClass>(k) == f.cls) ? p_own : (1.0 - p_own) / (kNumObjectClasses - 1);
      score[k] = p_radar + cfg_.camera_weight * ev[k] / ev_sum;
    }
    score[static_cast<std::size_t>(ObjectClass::kUnknown)] = -1.0;
    const auto best =
        static_cast<std::size_t>(std::max_element(score.begin(), score.end()) - score.begin());
    // Keep "static" only if the camera has never reported this track.
    if (f.cls != ObjectClass::kStatic || ev_sum > 1.5) {
      f.cls = static_cast<ObjectClass>(best);
      f.class_confidence = score[best] / (1.0 + cfg_.camera_weight);
    }
  }
  // Forget tracks that no longer exist.
  for (auto it = camera_evidence_.begin(); it != camera_evidence_.end();) {
    const bool alive = std::any_of(tracks.begin(), tracks.end(),
                                   [&](const Track& t) { return t.id == it->first; });
    it = alive ? std::next(it) : camera_evidence_.erase(it);
  }
  return out;
}

}  // namespace radar
