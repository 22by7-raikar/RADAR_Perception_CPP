#pragma once

#include <array>
#include <map>
#include <vector>

#include "radar/camera.hpp"
#include "radar/tracker.hpp"

namespace radar {

struct FusionConfig {
  double min_iou = 0.15;         ///< association gate on projected-track vs camera box
  bool compensate_time = true;   ///< predict tracks to the camera timestamp
  double radar_height_m = 0.6;   ///< radar above ground (projected boxes start at the ground)
  double camera_weight = 2.0;    ///< camera class evidence relative to radar
  double evidence_decay = 0.95;  ///< per camera frame, for the per-track camera evidence
};

/// One radar track, optionally confirmed by a camera detection.
struct FusedObject {
  int track_id = 0;
  int camera_index = -1;  ///< index into the camera detections, -1 if radar only
  ObjectClass cls = ObjectClass::kUnknown;
  double class_confidence = 0.0;
  Vec2 position_sensor = Vec2::Zero();  ///< at the camera timestamp
  Vec2 velocity_world = Vec2::Zero();
  BoundingBox projected;  ///< track box projected into the image
  double iou = 0.0;
};

/// Late fusion of radar tracks with camera 2D detections.
///
/// Each confirmed track is predicted to the camera timestamp (constant
/// velocity), moved into the radar frame at that time using the sensor pose
/// at the camera timestamp, given a 3D box (track extent and class prior),
/// projected through the calibrated extrinsic, and associated to camera boxes
/// by Hungarian assignment on 1 - IoU.
///
/// Classes are fused per track over time: camera detections accumulate decaying
/// class evidence, which is combined with the radar class probability. A single
/// misclassified camera frame cannot flip an established track, while a
/// pedestrian the radar called static is corrected after a few frames. Tracks
/// the radar labelled static are still projected, since that label is where
/// the radar is least reliable.
class RadarCameraFusion {
 public:
  RadarCameraFusion(const CameraModel& camera, const FusionConfig& cfg = {});

  std::vector<FusedObject> fuse(const std::vector<Track>& tracks, double radar_time_s,
                                const SensorState& sensor_at_camera_time, double camera_time_s,
                                const std::vector<CameraDetection>& detections);

  void set_camera(const CameraModel& camera) { camera_ = camera; }
  const FusionConfig& config() const { return cfg_; }

 private:
  CameraModel camera_;
  FusionConfig cfg_;
  std::map<int, std::array<double, kNumObjectClasses>> camera_evidence_;
};

}  // namespace radar
