#pragma once

#include <vector>

#include "radar/camera.hpp"

namespace radar {

/// A calibration target seen by both sensors: its 3D position from the radar
/// (x, y measured, z from the known target height) and its pixel in the image.
struct Correspondence {
  Vec3 p_radar = Vec3::Zero();
  Eigen::Vector2d uv = Eigen::Vector2d::Zero();
};

struct CalibrationResult {
  Pose3 T_cam_from_radar = Pose3::Identity();
  double initial_rms_px = 0.0;
  double final_rms_px = 0.0;
  int iterations = 0;
  bool converged = false;
};

/// Radar-to-camera extrinsic calibration from target correspondences.
///
/// Levenberg-Marquardt on SE(3) minimising reprojection error, with the
/// update applied as a left perturbation T <- exp(delta) T. Starting from a
/// rough mounting guess (CAD drawing, tape measure), it typically converges
/// in a handful of iterations.
CalibrationResult calibrate_extrinsic(const CameraIntrinsics& K,
                                      const std::vector<Correspondence>& data, const Pose3& initial,
                                      int max_iterations = 50);

/// Rotation angle (radians) and translation distance between two poses.
struct PoseError {
  double rotation_rad = 0.0;
  double translation_m = 0.0;
};
PoseError pose_difference(const Pose3& a, const Pose3& b);

}  // namespace radar
