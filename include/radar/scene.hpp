#pragma once

#include <string>
#include <vector>

#include "radar/geometry.hpp"
#include "radar/simulator.hpp"

namespace radar {

enum class ObjectClass : int { kUnknown = 0, kPedestrian = 1, kVehicle = 2, kStatic = 3 };
inline constexpr int kNumObjectClasses = 4;
const char* to_string(ObjectClass c);

/// A point reflector rigidly attached to an object.
struct Scatterer {
  Vec2 offset = Vec2::Zero();  ///< body frame: x along heading, y left
  Vec2 normal = Vec2::Zero();  ///< outward body-frame normal; zero = visible from all sides
  double rcs_dbsm = 0.0;
  /// Periodic motion along the heading (e.g. limbs): v += amp sin(2 pi f t + phase).
  double oscillation_amp_mps = 0.0;
  double oscillation_hz = 0.0;
  double oscillation_phase = 0.0;
};

/// Rigid object moving with constant speed and turn rate (CTRV).
struct SceneObject {
  int id = 0;
  ObjectClass cls = ObjectClass::kUnknown;
  Vec2 position = Vec2::Zero();  ///< world frame, object centre
  double heading = 0.0;
  double speed = 0.0;
  double yaw_rate = 0.0;
  double length = 0.0;
  double width = 0.0;
  std::vector<Scatterer> scatterers;

  Vec2 velocity() const { return speed * Vec2(std::cos(heading), std::sin(heading)); }

  static SceneObject vehicle(int id, Vec2 position, double heading, double speed,
                             double yaw_rate = 0.0, double length = 4.6, double width = 1.9);
  static SceneObject pedestrian(int id, Vec2 position, double heading, double speed);
  static SceneObject pole(int id, Vec2 position, double rcs_dbsm = 8.0);
  /// Line of stationary reflectors (guard rail, fence, wall) from `a` to `b`.
  static SceneObject barrier(int id, Vec2 a, Vec2 b, double spacing_m = 4.0, double rcs_dbsm = 3.0);
};

/// Ego platform state in the world frame. Constant acceleration and turn rate.
struct EgoState {
  Pose2D pose;
  double speed = 0.0;  ///< along heading
  double yaw_rate = 0.0;
  double accel = 0.0;

  Vec2 velocity() const { return speed * Vec2(std::cos(pose.yaw), std::sin(pose.yaw)); }
};

struct GroundTruthObject {
  int id = 0;
  ObjectClass cls = ObjectClass::kUnknown;
  Vec2 position_world = Vec2::Zero();
  Vec2 velocity_world = Vec2::Zero();
  Vec2 position_sensor = Vec2::Zero();
  double heading = 0.0;
  double length = 0.0;
  double width = 0.0;
  std::size_t num_scatterers = 0;  ///< scatterers that produced an echo this frame
};

struct SceneConfig {
  Pose2D sensor_mount;                            ///< sensor pose in the ego body frame
  double field_of_view_rad = 75.0 * kPi / 180.0;  ///< half-angle
  double min_range_m = 1.0;
  double max_range_m = 250.0;
};

struct SceneFrame {
  double time_s = 0.0;
  EgoState ego;
  Pose2D sensor_pose;                      ///< world frame
  Vec2 sensor_velocity = Vec2::Zero();     ///< world frame
  std::vector<Target> targets;             ///< echoes in sensor coordinates, input to Simulator
  std::vector<GroundTruthObject> objects;  ///< objects with at least one echo
};

/// World model that turns objects and ego motion into per-frame radar targets
/// (range, radial velocity, azimuth) and ground truth.
class Scene {
 public:
  explicit Scene(const SceneConfig& cfg = {}, const EgoState& ego = {});

  void add(SceneObject obj) { objects_.push_back(std::move(obj)); }
  void step(double dt);
  SceneFrame snapshot() const;

  double time() const { return time_; }
  const EgoState& ego() const { return ego_; }
  EgoState& ego() { return ego_; }
  const std::vector<SceneObject>& objects() const { return objects_; }
  const SceneConfig& config() const { return cfg_; }

 private:
  SceneConfig cfg_;
  EgoState ego_;
  std::vector<SceneObject> objects_;
  double time_ = 0.0;
};

/// Constant turn rate and velocity propagation of a planar pose.
void propagate_ctrv(Vec2& position, double& heading, double speed, double yaw_rate, double dt);

}  // namespace radar
