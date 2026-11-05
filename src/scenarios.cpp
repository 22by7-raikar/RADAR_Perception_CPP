#include "radar/scenarios.hpp"

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

Scene ego_drive(const SceneConfig& cfg) {
  EgoState ego;
  ego.speed = 8.0;
  ego.yaw_rate = 0.03;
  ego.accel = 0.3;
  Scene s(cfg, ego);
  int id = 1;
  s.add(SceneObject::vehicle(id++, {60.0, 3.5}, 0.0, 11.0));             // lead vehicle
  s.add(SceneObject::vehicle(id++, {160.0, -3.5}, 180.0 * kDeg, 10.0));  // oncoming
  s.add(SceneObject::vehicle(id++, {90.0, -40.0}, 90.0 * kDeg, 7.0));    // crossing
  s.add(SceneObject::pedestrian(id++, {45.0, 8.0}, -100.0 * kDeg, 1.3));
  add_environment(s, id);
  for (double x : {140.0, 170.0, 200.0}) s.add(SceneObject::pole(id++, {x, 12.0}));
  s.add(SceneObject::barrier(id++, {10.0, 14.0}, {220.0, 20.0}, 6.0));
  return s;
}

}  // namespace radar::scenarios
