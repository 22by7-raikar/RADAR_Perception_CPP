#pragma once

#include <Eigen/Core>
#include <array>
#include <vector>

#include "radar/clustering.hpp"
#include "radar/geometry.hpp"
#include "radar/scene.hpp"

namespace radar {

/// Sensor pose and velocity in the tracking (world / odometry) frame.
struct SensorState {
  Pose2D pose;
  Vec2 velocity = Vec2::Zero();
};

struct TrackerConfig {
  double accel_std_mps2 = 2.5;  ///< white-acceleration process noise
  double range_std_m = 0.4;
  double azimuth_std_rad = 1.5 * kPi / 180.0;
  double extent_std_per_m = 0.35;  ///< extra centroid std per metre of cluster extent
  double radial_velocity_std_mps = 0.4;
  double gate_chi2 = 13.8;  ///< 3-DoF gate, ~99.7 %
  int confirm_hits = 3;
  int max_tentative_misses = 1;
  int max_confirmed_misses = 6;
  double init_los_velocity_std_mps = 1.0;
  double init_tangential_velocity_std_mps = 8.0;
  /// Gate for "part of the same extended object": along the track heading,
  /// across it, and in radial velocity. Used to absorb leftover clusters and
  /// to merge duplicate tracks.
  double absorb_longitudinal_m = 5.0;
  double absorb_lateral_m = 1.6;
  double absorb_velocity_mps = 1.0;
  double merge_velocity_mps = 1.5;
  double merge_velocity_chi2 = 9.2;  ///< 2-DoF, 99 %: velocities statistically equal
};

struct ClassifierConfig {
  double static_speed_mps = 0.6;
  double static_velocity_chi2 = 6.0;  ///< 2-DoF, 95 %: young track velocity consistent with zero
  double static_max_speed_mps = 2.0;
  /// After this long, "static" is judged from net displacement over the track
  /// life, which angle noise on distant objects does not accumulate.
  double static_min_history_s = 1.0;
  double static_mean_speed_mps = 0.4;
  /// Static objects have no Doppler spread and zero ego-compensated radial velocity.
  double static_max_spread_mps = 0.5;
  double static_max_radial_mps = 0.5;
  double vehicle_rcs_dbsm = 8.0;
  double vehicle_extent_m = 1.8;
  double pedestrian_rcs_dbsm = 5.0;
  double pedestrian_spread_mps = 0.9;
  double pedestrian_max_speed_mps = 3.5;
  double feature_smoothing = 0.3;  ///< EMA weight of the newest observation
  double evidence_decay = 0.9;
  double decision_threshold = 1.5;
};

enum class TrackStatus { kTentative, kConfirmed };

struct Track {
  using State = Eigen::Vector4d;  ///< [px, py, vx, vy], world frame
  using Covariance = Eigen::Matrix4d;

  int id = 0;
  TrackStatus status = TrackStatus::kTentative;
  State x = State::Zero();
  Covariance P = Covariance::Identity();
  int hits = 0;
  int consecutive_misses = 0;
  int age = 0;
  double last_update_s = 0.0;
  double first_update_s = 0.0;
  Eigen::Vector2d first_position = Eigen::Vector2d::Zero();
  double first_position_std_m = 0.0;

  // Classification
  ObjectClass cls = ObjectClass::kUnknown;
  double class_confidence = 0.0;
  std::array<double, kNumObjectClasses> evidence{};
  double extent_m = 0.0;
  double rcs_dbsm = 0.0;
  double velocity_spread_mps = 0.0;
  double radial_world_mps = 0.0;  ///< |radial velocity| relative to the static world (EMA)
  bool has_features = false;

  Vec2 position() const { return x.head<2>(); }
  Vec2 velocity() const { return x.tail<2>(); }
  bool confirmed() const { return status == TrackStatus::kConfirmed; }
};

/// Multi-object tracker: constant-velocity EKF per track, measurements of
/// cluster position plus radial velocity, global-nearest-neighbour association
/// (Hungarian on gated Mahalanobis cost), M-of-N confirmation and miss-based
/// deletion, and per-track class evidence accumulation.
///
/// Tracks live in a world/odometry frame; the sensor pose and velocity at each
/// update are inputs, so ego motion is compensated exactly.
class MultiTargetTracker {
 public:
  explicit MultiTargetTracker(const TrackerConfig& cfg = {}, const ClassifierConfig& cls = {});

  void update(double time_s, const std::vector<Cluster>& clusters, const SensorState& sensor);

  const std::vector<Track>& tracks() const { return tracks_; }
  std::vector<Track> confirmed_tracks() const;
  void reset();

 private:
  struct Measurement {
    Eigen::Vector3d z;
    Eigen::Matrix3d R;
  };

  Measurement make_measurement(const Cluster& c, const SensorState& s) const;
  void predict(Track& t, double dt) const;
  /// Predicted measurement and its Jacobian.
  Eigen::Vector3d h(const Track& t, const SensorState& s, Eigen::Matrix<double, 3, 4>* H) const;
  void correct(Track& t, const Measurement& m, const SensorState& s) const;
  Track initiate(const Measurement& m, const SensorState& s, double time_s);
  void classify(Track& t, const Cluster& c, const SensorState& sensor, double time_s) const;
  void merge_duplicates();
  bool within_object_extent(const Track& t, const Vec2& p) const;

  TrackerConfig cfg_;
  ClassifierConfig cls_cfg_;
  std::vector<Track> tracks_;
  int next_id_ = 1;
  double last_time_s_ = 0.0;
  bool initialised_ = false;
};

}  // namespace radar
