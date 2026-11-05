#include "radar/scene.hpp"

#include <cmath>

namespace radar {

const char* to_string(ObjectClass c) {
  switch (c) {
    case ObjectClass::kPedestrian:
      return "pedestrian";
    case ObjectClass::kVehicle:
      return "vehicle";
    case ObjectClass::kStatic:
      return "static";
    case ObjectClass::kUnknown:
      break;
  }
  return "unknown";
}

SceneObject SceneObject::vehicle(int id, Vec2 position, double heading, double speed,
                                 double yaw_rate, double length, double width) {
  SceneObject o;
  o.id = id;
  o.cls = ObjectClass::kVehicle;
  o.position = position;
  o.heading = heading;
  o.speed = speed;
  o.yaw_rate = yaw_rate;
  o.length = length;
  o.width = width;
  const double hl = 0.5 * length, hw = 0.5 * width;
  const double r = std::sqrt(0.5);
  // Corners are strong (dihedral-like), side and end points weaker.
  o.scatterers = {
      {{hl, hw}, {r, r}, 12.0},        {{hl, -hw}, {r, -r}, 12.0},
      {{-hl, hw}, {-r, r}, 12.0},      {{-hl, -hw}, {-r, -r}, 12.0},
      {{hl, 0.0}, {1, 0}, 10.0},       {{-hl, 0.0}, {-1, 0}, 10.0},
      {{0.3 * hl, hw}, {0, 1}, 6.0},   {{-0.5 * hl, hw}, {0, 1}, 6.0},
      {{0.3 * hl, -hw}, {0, -1}, 6.0}, {{-0.5 * hl, -hw}, {0, -1}, 6.0},
  };
  return o;
}

SceneObject SceneObject::pedestrian(int id, Vec2 position, double heading, double speed) {
  SceneObject o;
  o.id = id;
  o.cls = ObjectClass::kPedestrian;
  o.position = position;
  o.heading = heading;
  o.speed = speed;
  o.length = 0.4;
  o.width = 0.5;
  const double gait_hz = 0.9 + 0.4 * speed;  // stride frequency grows with speed
  const double swing = 0.6 * speed + 0.8;    // limb speed relative to torso
  o.scatterers = {
      {{0.0, 0.0}, Vec2::Zero(), -1.0},                         // torso
      {{0.1, 0.1}, Vec2::Zero(), -8.0, swing, gait_hz, 0.0},    // limb
      {{-0.1, -0.1}, Vec2::Zero(), -8.0, swing, gait_hz, kPi},  // opposite limb
  };
  return o;
}

SceneObject SceneObject::pole(int id, Vec2 position, double rcs_dbsm) {
  SceneObject o;
  o.id = id;
  o.cls = ObjectClass::kStatic;
  o.position = position;
  o.length = o.width = 0.3;
  o.scatterers = {{Vec2::Zero(), Vec2::Zero(), rcs_dbsm}};
  return o;
}

SceneObject SceneObject::barrier(int id, Vec2 a, Vec2 b, double spacing_m, double rcs_dbsm) {
  SceneObject o;
  o.id = id;
  o.cls = ObjectClass::kStatic;
  o.position = 0.5 * (a + b);
  const Vec2 d = b - a;
  o.heading = std::atan2(d.y(), d.x());
  o.length = d.norm();
  o.width = 0.2;
  const auto n = static_cast<int>(std::floor(o.length / spacing_m));
  for (int i = 0; i <= n; ++i) {
    const double s = -0.5 * o.length + spacing_m * i;
    o.scatterers.push_back({{s, 0.0}, Vec2::Zero(), rcs_dbsm});
  }
  return o;
}

void propagate_ctrv(Vec2& position, double& heading, double speed, double yaw_rate, double dt) {
  if (std::abs(yaw_rate) < 1e-9) {
    position += speed * dt * Vec2(std::cos(heading), std::sin(heading));
  } else {
    const double h1 = heading + yaw_rate * dt;
    position += (speed / yaw_rate) *
                Vec2(std::sin(h1) - std::sin(heading), std::cos(heading) - std::cos(h1));
    heading = h1;
  }
  heading = wrap_angle(heading);
}

Scene::Scene(const SceneConfig& cfg, const EgoState& ego) : cfg_(cfg), ego_(ego) {}

void Scene::step(double dt) {
  for (SceneObject& o : objects_) propagate_ctrv(o.position, o.heading, o.speed, o.yaw_rate, dt);
  // Ego: integrate with the mid-interval speed so constant acceleration is exact on straights.
  Vec2 p = ego_.pose.translation();
  double yaw = ego_.pose.yaw;
  const double v_mid = ego_.speed + 0.5 * ego_.accel * dt;
  propagate_ctrv(p, yaw, v_mid, ego_.yaw_rate, dt);
  ego_.pose = {p.x(), p.y(), yaw};
  ego_.speed += ego_.accel * dt;
  time_ += dt;
}

SceneFrame Scene::snapshot() const {
  SceneFrame f;
  f.time_s = time_;
  f.ego = ego_;
  f.sensor_pose = ego_.pose.compose(cfg_.sensor_mount);
  const Vec2 mount_world = ego_.pose.rotation_matrix() * cfg_.sensor_mount.translation();
  f.sensor_velocity = ego_.velocity() + cross_z(ego_.yaw_rate, mount_world);
  const Vec2 sensor_pos = f.sensor_pose.translation();

  for (const SceneObject& o : objects_) {
    const Mat2 rot = rotation(o.heading);
    const Vec2 v_body = o.velocity();
    GroundTruthObject gt{
        o.id,      o.cls,    o.position, v_body, f.sensor_pose.to_local(o.position),
        o.heading, o.length, o.width,    0};

    for (const Scatterer& s : o.scatterers) {
      const Vec2 p = o.position + rot * s.offset;
      const Vec2 to_sensor = sensor_pos - p;
      if (s.normal.squaredNorm() > 0.0 && (rot * s.normal).dot(to_sensor) <= 0.0) continue;

      const Vec2 local = f.sensor_pose.to_local(p);
      const double range = local.norm();
      const double az = std::atan2(local.y(), local.x());
      if (range < cfg_.min_range_m || range > cfg_.max_range_m) continue;
      if (std::abs(az) > cfg_.field_of_view_rad) continue;

      Vec2 v = v_body + cross_z(o.yaw_rate, rot * s.offset);
      if (s.oscillation_amp_mps != 0.0) {
        const double osc = s.oscillation_amp_mps *
                           std::sin(2.0 * kPi * s.oscillation_hz * time_ + s.oscillation_phase);
        v += osc * Vec2(std::cos(o.heading), std::sin(o.heading));
      }
      // Radial velocity is rotation invariant: d . (v_point - v_sensor) / |d|.
      const Vec2 d = p - sensor_pos;
      const double radial = d.dot(v - f.sensor_velocity) / d.norm();
      f.targets.push_back({range, radial, az, s.rcs_dbsm});
      ++gt.num_scatterers;
    }
    if (gt.num_scatterers > 0) f.objects.push_back(gt);
  }
  return f;
}

}  // namespace radar
