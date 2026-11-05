#include "radar/scenarios.hpp"

#include <cmath>
#include <random>
#include <utility>
#include <vector>

namespace radar::scenarios {

namespace {

constexpr double kDeg = kPi / 180.0;

void add_environment(Scene& s, int& id) {
  for (double x : {20.0, 50.0, 80.0, 110.0}) s.add(SceneObject::pole(id++, {x, -7.0}));
  for (double x : {35.0, 95.0}) s.add(SceneObject::pole(id++, {x, 9.0}));
  s.add(SceneObject::barrier(id++, {10.0, -10.0}, {170.0, -10.0}, 5.0));
}

}  // namespace

Scene intersection(const SceneConfig& cfg) {
  Scene s(cfg);
  int id = 1;
  s.add(SceneObject::vehicle(id++, {70.0, -30.0}, 90.0 * kDeg, 8.0));        // crossing left
  s.add(SceneObject::vehicle(id++, {130.0, 3.5}, 180.0 * kDeg, 12.0));       // oncoming
  s.add(SceneObject::vehicle(id++, {25.0, -3.5}, 0.0, 6.0));                 // receding
  s.add(SceneObject::vehicle(id++, {55.0, 25.0}, -90.0 * kDeg, 7.0, 0.12));  // turning
  s.add(SceneObject::pedestrian(id++, {22.0, 9.0}, -90.0 * kDeg, 1.4));      // crossing
  s.add(SceneObject::pedestrian(id++, {40.0, -7.0}, 35.0 * kDeg, 1.1));
  add_environment(s, id);
  return s;
}

Scene ego_drive(const SceneConfig& cfg, double duration_s) {
  EgoState ego;
  ego.speed = 8.0;
  ego.yaw_rate = 0.03;
  ego.accel = 0.3;
  SceneConfig c = cfg;
  if (c.sensor_mount.translation().norm() == 0.0) c.sensor_mount = {3.6, 0.0, 0.0};  // bumper
  Scene s(c, ego);

  // Lay out roadside poles along the ego's future path so that stationary
  // reflectors stay in view for the whole drive.
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> rcs(2.0, 10.0), lateral(8.0, 11.0);
  EgoState probe = ego;
  Vec2 p = probe.pose.translation();
  double yaw = probe.pose.yaw, travelled = 0.0, next = 15.0;
  const double dt = 0.05;
  int id = 100;
  std::vector<std::pair<Vec2, double>> path;  // (point, heading) every metre or so
  for (double t = 0.0; t < duration_s + 15.0; t += dt) {
    const Vec2 before = p;
    const double v = probe.speed + 0.5 * probe.accel * dt;
    propagate_ctrv(p, yaw, v, probe.yaw_rate, dt);
    probe.speed += probe.accel * dt;
    travelled += (p - before).norm();
    path.emplace_back(p, yaw);
    if (travelled >= next) {
      const Vec2 left(-std::sin(yaw), std::cos(yaw));
      s.add(SceneObject::pole(id++, p + lateral(rng) * left, rcs(rng)));
      s.add(SceneObject::pole(id++, p - lateral(rng) * left, rcs(rng)));
      next += 12.0;
    }
  }
  auto at = [&](double metres) {
    // Path point roughly `metres` along the route.
    double acc = 0.0;
    for (std::size_t i = 1; i < path.size(); ++i) {
      acc += (path[i].first - path[i - 1].first).norm();
      if (acc >= metres) return path[i];
    }
    return path.back();
  };

  int oid = 1;
  {
    const auto [q, h] = at(45.0);  // lead vehicle in the ego lane
    s.add(SceneObject::vehicle(oid++, q, h, 11.0, 0.03));
  }
  {
    const auto [q, h] = at(230.0);  // oncoming, adjacent lane
    const Vec2 left(-std::sin(h), std::cos(h));
    s.add(SceneObject::vehicle(oid++, q + 3.5 * left, wrap_angle(h + kPi), 12.0, -0.035));
  }
  {
    const auto [q, h] = at(140.0);  // crossing traffic
    const Vec2 left(-std::sin(h), std::cos(h));
    s.add(SceneObject::vehicle(oid++, q - 35.0 * left, wrap_angle(h + kPi / 2.0), 7.0));
  }
  {
    const auto [q, h] = at(90.0);  // pedestrian crossing the road
    const Vec2 left(-std::sin(h), std::cos(h));
    s.add(SceneObject::pedestrian(oid++, q + 9.0 * left, wrap_angle(h - kPi / 2.0), 1.3));
  }
  return s;
}

}  // namespace radar::scenarios
