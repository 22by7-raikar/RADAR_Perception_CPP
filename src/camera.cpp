#include "radar/camera.hpp"

#include <algorithm>
#include <cmath>

namespace radar {

Eigen::Matrix3d radar_to_camera_axes() {
  Eigen::Matrix3d r;
  // camera x (right) = -radar y, camera y (down) = -radar z, camera z = radar x
  r << 0, -1, 0, 0, 0, -1, 1, 0, 0;
  return r;
}

Pose3 make_camera_extrinsic(const Vec3& position, double roll, double pitch, double yaw) {
  // Misalignment applied in the radar frame: R_mis = Rz(yaw) Ry(pitch) Rx(roll).
  const Eigen::Matrix3d mis =
      (Eigen::AngleAxisd(yaw, Vec3::UnitZ()) * Eigen::AngleAxisd(pitch, Vec3::UnitY()) *
       Eigen::AngleAxisd(roll, Vec3::UnitX()))
          .toRotationMatrix();
  // Camera pose in the radar frame, then invert to get T_cam_from_radar.
  Pose3 radar_from_cam = Pose3::Identity();
  radar_from_cam.linear() = mis * radar_to_camera_axes().transpose();
  radar_from_cam.translation() = position;
  return radar_from_cam.inverse();
}

double BoundingBox::iou(const BoundingBox& o) const {
  const BoundingBox in{std::max(u0, o.u0), std::max(v0, o.v0), std::min(u1, o.u1),
                       std::min(v1, o.v1)};
  const double i = in.area();
  const double u = area() + o.area() - i;
  return u > 0.0 ? i / u : 0.0;
}

bool CameraModel::project(const Vec3& p_radar, Eigen::Vector2d& uv) const {
  const Vec3 p = T_cam_from_radar * p_radar;
  if (p.z() < 0.5) return false;
  uv << K.fx * p.x() / p.z() + K.cx, K.fy * p.y() / p.z() + K.cy;
  return true;
}

bool CameraModel::project_box(const Vec2& centre, double heading, double length, double width,
                              double height, double ground_z, BoundingBox& box) const {
  const double c = std::cos(heading), s = std::sin(heading);
  box = {1e18, 1e18, -1e18, -1e18};
  for (int i = 0; i < 8; ++i) {
    const double lx = (i & 1 ? 0.5 : -0.5) * length;
    const double ly = (i & 2 ? 0.5 : -0.5) * width;
    const double z = ground_z + (i & 4 ? height : 0.0);
    const Vec3 p(centre.x() + c * lx - s * ly, centre.y() + s * lx + c * ly, z);
    Eigen::Vector2d uv;
    if (!project(p, uv)) return false;
    box.u0 = std::min(box.u0, uv.x());
    box.v0 = std::min(box.v0, uv.y());
    box.u1 = std::max(box.u1, uv.x());
    box.v1 = std::max(box.v1, uv.y());
  }
  // Clip to the image; reject boxes entirely outside.
  box.u0 = std::max(box.u0, 0.0);
  box.v0 = std::max(box.v0, 0.0);
  box.u1 = std::min(box.u1, static_cast<double>(K.width));
  box.v1 = std::min(box.v1, static_cast<double>(K.height));
  return box.u1 > box.u0 && box.v1 > box.v0;
}

ObjectSize nominal_size(ObjectClass c) {
  switch (c) {
    case ObjectClass::kPedestrian:
      return {0.5, 0.6, 1.75};
    case ObjectClass::kVehicle:
      return {4.6, 1.9, 1.5};
    case ObjectClass::kStatic:
      return {0.3, 0.3, 3.0};
    case ObjectClass::kUnknown:
      break;
  }
  return {1.5, 1.5, 1.5};
}

CameraSimulator::CameraSimulator(const CameraModel& camera, double radar_height_m,
                                 const CameraDetectorConfig& cfg)
    : camera_(camera), radar_height_(radar_height_m), cfg_(cfg), rng_(cfg.seed) {}

std::vector<CameraDetection> CameraSimulator::detect(const SceneFrame& frame) {
  std::vector<CameraDetection> out;
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  std::normal_distribution<double> px(0.0, cfg_.pixel_noise);

  for (const GroundTruthObject& g : frame.objects) {
    if (g.cls != ObjectClass::kVehicle && g.cls != ObjectClass::kPedestrian) continue;
    const ObjectSize sz = nominal_size(g.cls);
    const double heading = wrap_angle(g.heading - frame.sensor_pose.yaw);
    BoundingBox b;
    if (!camera_.project_box(g.position_sensor, heading, g.length, g.width, sz.height,
                             -radar_height_, b)) {
      continue;
    }
    if (b.v1 - b.v0 < cfg_.min_box_height_px || u01(rng_) < cfg_.miss_probability) continue;
    CameraDetection d;
    d.box = {b.u0 + px(rng_), b.v0 + px(rng_), b.u1 + px(rng_), b.v1 + px(rng_)};
    d.truth_id = g.id;
    d.cls = g.cls;
    if (u01(rng_) < cfg_.class_confusion) {
      d.cls = g.cls == ObjectClass::kVehicle ? ObjectClass::kPedestrian : ObjectClass::kVehicle;
    }
    d.score = 0.6 + 0.4 * u01(rng_);
    out.push_back(d);
  }

  // Poisson-distributed false positives with random size and position.
  std::poisson_distribution<int> nfp(cfg_.false_positives_per_frame);
  const int fps = nfp(rng_);
  for (int i = 0; i < fps; ++i) {
    const double w = 30.0 + 200.0 * u01(rng_), h = 30.0 + 150.0 * u01(rng_);
    const double u = u01(rng_) * (camera_.K.width - w), v = 300.0 + u01(rng_) * 400.0;
    CameraDetection d;
    d.box = {u, v, u + w, v + h};
    d.cls = u01(rng_) < 0.5 ? ObjectClass::kVehicle : ObjectClass::kPedestrian;
    d.score = 0.3 + 0.4 * u01(rng_);
    out.push_back(d);
  }
  return out;
}

}  // namespace radar
