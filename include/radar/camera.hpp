#pragma once

#include <Eigen/Geometry>
#include <cstdint>
#include <random>
#include <vector>

#include "radar/scene.hpp"

namespace radar {

using Vec3 = Eigen::Vector3d;
using Pose3 = Eigen::Isometry3d;

/// Rotation from radar axes (x forward, y left, z up) to camera axes
/// (x right, y down, z forward), before any mounting misalignment.
Eigen::Matrix3d radar_to_camera_axes();

/// T_cam_from_radar for a camera at `position` (radar frame) with small
/// roll/pitch/yaw misalignments (radians) relative to the radar axes.
Pose3 make_camera_extrinsic(const Vec3& position, double roll, double pitch, double yaw);

struct CameraIntrinsics {
  double fx = 1000.0;
  double fy = 1000.0;
  double cx = 960.0;
  double cy = 540.0;
  int width = 1920;
  int height = 1080;
};

struct BoundingBox {
  double u0 = 0.0, v0 = 0.0, u1 = 0.0, v1 = 0.0;
  double area() const { return std::max(0.0, u1 - u0) * std::max(0.0, v1 - v0); }
  Eigen::Vector2d center() const { return {0.5 * (u0 + u1), 0.5 * (v0 + v1)}; }
  double iou(const BoundingBox& o) const;
};

/// Pinhole camera with its extrinsic relative to the radar.
struct CameraModel {
  CameraIntrinsics K;
  Pose3 T_cam_from_radar = Pose3::Identity();

  /// Projects a radar-frame point. False if it is behind the camera.
  bool project(const Vec3& p_radar, Eigen::Vector2d& uv) const;
  /// Image box of an upright 3D box (radar frame: centre on the ground at
  /// z = ground_z, heading in the radar frame). False if not visible.
  bool project_box(const Vec2& centre, double heading, double length, double width, double height,
                   double ground_z, BoundingBox& box) const;
};

/// Typical object dimensions for the camera detector and for projecting tracks.
struct ObjectSize {
  double length, width, height;
};
ObjectSize nominal_size(ObjectClass c);

struct CameraDetection {
  BoundingBox box;
  ObjectClass cls = ObjectClass::kUnknown;
  double score = 0.0;
  int truth_id = -1;  ///< ground-truth object (-1 for a false positive); for evaluation only
};

/// 2D object detector on a camera image: projects ground-truth boxes, adds
/// corner noise, misses, false positives and class confusion.
struct CameraDetectorConfig {
  double pixel_noise = 4.0;
  double miss_probability = 0.05;
  double false_positives_per_frame = 0.3;
  double class_confusion = 0.03;
  double min_box_height_px = 12.0;  ///< smaller objects are not detected
  std::uint32_t seed = 21;
};

class CameraSimulator {
 public:
  CameraSimulator(const CameraModel& camera, double radar_height_m,
                  const CameraDetectorConfig& cfg = {});

  /// Detections for objects in `frame` (positions in the sensor frame of the
  /// frame's snapshot). Only pedestrians and vehicles are detected.
  std::vector<CameraDetection> detect(const SceneFrame& frame);

  const CameraModel& camera() const { return camera_; }

 private:
  CameraModel camera_;
  double radar_height_;
  CameraDetectorConfig cfg_;
  std::mt19937 rng_;
};

}  // namespace radar
