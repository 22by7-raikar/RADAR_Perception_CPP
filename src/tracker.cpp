#include "radar/tracker.hpp"

#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>

#include "radar/assignment.hpp"

namespace radar {

namespace {
constexpr double kForbidden = 1e9;
}

MultiTargetTracker::MultiTargetTracker(const TrackerConfig& cfg, const ClassifierConfig& cls)
    : cfg_(cfg), cls_cfg_(cls) {}

void MultiTargetTracker::reset() {
  tracks_.clear();
  next_id_ = 1;
  initialised_ = false;
}

std::vector<Track> MultiTargetTracker::confirmed_tracks() const {
  std::vector<Track> out;
  for (const Track& t : tracks_) {
    if (t.confirmed()) out.push_back(t);
  }
  return out;
}

MultiTargetTracker::Measurement MultiTargetTracker::make_measurement(const Cluster& c,
                                                                     const SensorState& s) const {
  // Polar measurement noise mapped to Cartesian (sensor frame), inflated by the
  // cluster extent because the centroid of an extended object wanders between
  // frames as different scatterers are resolved.
  const double r = std::max(0.1, c.centroid.norm());
  const double az = std::atan2(c.centroid.y(), c.centroid.x());
  Mat2 J;
  J << std::cos(az), -r * std::sin(az), std::sin(az), r * std::cos(az);
  Mat2 polar = Mat2::Zero();
  polar(0, 0) = cfg_.range_std_m * cfg_.range_std_m;
  polar(1, 1) = cfg_.azimuth_std_rad * cfg_.azimuth_std_rad;
  const double ext = cfg_.extent_std_per_m * c.extent_major;
  const Mat2 cov_sensor = J * polar * J.transpose() + ext * ext * Mat2::Identity();
  const Mat2 rot = s.pose.rotation_matrix();

  Measurement m;
  m.z.head<2>() = s.pose.to_parent(c.centroid);
  m.z(2) = c.radial_velocity;
  m.R.setZero();
  m.R.topLeftCorner<2, 2>() = rot * cov_sensor * rot.transpose();
  const double sv = cfg_.radial_velocity_std_mps + 0.5 * c.velocity_spread;
  m.R(2, 2) = sv * sv;
  return m;
}

void MultiTargetTracker::predict(Track& t, double dt) const {
  if (dt <= 0.0) return;
  Eigen::Matrix4d F = Eigen::Matrix4d::Identity();
  F(0, 2) = F(1, 3) = dt;
  const double q = cfg_.accel_std_mps2 * cfg_.accel_std_mps2;
  const double dt2 = dt * dt, dt3 = dt2 * dt, dt4 = dt3 * dt;
  Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
  Q(0, 0) = Q(1, 1) = q * dt4 / 4.0;
  Q(0, 2) = Q(2, 0) = Q(1, 3) = Q(3, 1) = q * dt3 / 2.0;
  Q(2, 2) = Q(3, 3) = q * dt2;
  t.x = F * t.x;
  t.P = F * t.P * F.transpose() + Q;
}

Eigen::Vector3d MultiTargetTracker::h(const Track& t, const SensorState& s,
                                      Eigen::Matrix<double, 3, 4>* H) const {
  // Radial velocity d . (v - v_s) / |d| does not depend on sensor rotation.
  const Vec2 d = t.position() - s.pose.translation();
  const Vec2 u = t.velocity() - s.velocity;
  const double r = std::max(0.1, d.norm());
  Eigen::Vector3d z;
  z << t.x(0), t.x(1), d.dot(u) / r;
  if (H) {
    H->setZero();
    (*H)(0, 0) = (*H)(1, 1) = 1.0;
    const Vec2 dvr_dp = u / r - d * (d.dot(u) / (r * r * r));
    (*H)(2, 0) = dvr_dp.x();
    (*H)(2, 1) = dvr_dp.y();
    (*H)(2, 2) = d.x() / r;
    (*H)(2, 3) = d.y() / r;
  }
  return z;
}

void MultiTargetTracker::correct(Track& t, const Measurement& m, const SensorState& s) const {
  Eigen::Matrix<double, 3, 4> H;
  const Eigen::Vector3d y = m.z - h(t, s, &H);
  const Eigen::Matrix3d S = H * t.P * H.transpose() + m.R;
  const Eigen::Matrix<double, 4, 3> K = t.P * H.transpose() * S.inverse();
  t.x += K * y;
  // Joseph form keeps P symmetric positive definite.
  const Eigen::Matrix4d I_KH = Eigen::Matrix4d::Identity() - K * H;
  t.P = I_KH * t.P * I_KH.transpose() + K * m.R * K.transpose();
  t.P = 0.5 * (t.P + t.P.transpose());
}

Track MultiTargetTracker::initiate(const Measurement& m, const SensorState& s, double time_s) {
  Track t;
  t.id = next_id_++;
  t.last_update_s = time_s;
  t.first_update_s = time_s;
  t.hits = 1;
  t.age = 1;
  const Vec2 p = m.z.head<2>();
  const Vec2 d = p - s.pose.translation();
  const Vec2 u = d.normalized();
  const Vec2 tangent(-u.y(), u.x());
  // Only the line-of-sight velocity is observed: u . v = v_r + u . v_s.
  const Vec2 v = u * (m.z(2) + u.dot(s.velocity));
  t.x << p.x(), p.y(), v.x(), v.y();
  t.first_position = p;
  t.first_position_std_m = std::sqrt(m.R.topLeftCorner<2, 2>().trace());
  t.P.setZero();
  t.P.topLeftCorner<2, 2>() = m.R.topLeftCorner<2, 2>();
  const double sl = cfg_.init_los_velocity_std_mps, st = cfg_.init_tangential_velocity_std_mps;
  t.P.bottomRightCorner<2, 2>() =
      sl * sl * u * u.transpose() + st * st * tangent * tangent.transpose();
  return t;
}

void MultiTargetTracker::classify(Track& t, const Cluster& c, const SensorState& sensor,
                                  double time_s) const {
  const double a = cls_cfg_.feature_smoothing;
  // Radial velocity relative to the static world: zero for anything stationary.
  const Vec2 los = (sensor.pose.to_parent(c.centroid) - sensor.pose.translation()).normalized();
  const double vr_world = c.radial_velocity + los.dot(sensor.velocity);
  if (!t.has_features) {
    t.radial_world_mps = std::abs(vr_world);
    t.extent_m = c.extent_major;
    t.rcs_dbsm = c.rcs_dbsm;
    t.velocity_spread_mps = c.velocity_spread;
    t.has_features = true;
  } else {
    t.radial_world_mps = (1 - a) * t.radial_world_mps + a * std::abs(vr_world);
    t.extent_m = (1 - a) * t.extent_m + a * c.extent_major;
    t.rcs_dbsm = (1 - a) * t.rcs_dbsm + a * c.rcs_dbsm;
    t.velocity_spread_mps = (1 - a) * t.velocity_spread_mps + a * c.velocity_spread;
  }

  const double speed = t.velocity().norm();
  // Static if slow, or if the velocity estimate is statistically consistent
  // with zero (far tracks have noisy tangential velocity from angle noise).
  const Vec2 v = t.velocity();
  const Mat2 Pv = t.P.bottomRightCorner<2, 2>();
  const double v_chi2 = v.dot(Pv.ldlt().solve(v));
  const double history = time_s - t.first_update_s;
  bool looks_static;
  if (history >= cls_cfg_.static_min_history_s) {
    // Displacement beyond what the position uncertainty at either end explains.
    const double tol =
        2.0 * (t.first_position_std_m + std::sqrt(t.P.topLeftCorner<2, 2>().trace()));
    const double moved = std::max(0.0, (t.position() - t.first_position).norm() - tol);
    const double mean_speed = moved / history;
    looks_static = mean_speed < cls_cfg_.static_mean_speed_mps;
  } else {
    looks_static = speed < cls_cfg_.static_speed_mps || (v_chi2 < cls_cfg_.static_velocity_chi2 &&
                                                         speed < cls_cfg_.static_max_speed_mps);
  }
  // Doppler evidence of motion overrides kinematic "static": limbs spread the
  // Doppler of a walking person even when it moves across the line of sight.
  if (t.velocity_spread_mps > cls_cfg_.static_max_spread_mps ||
      t.radial_world_mps > cls_cfg_.static_max_radial_mps) {
    looks_static = false;
  }
  std::array<double, kNumObjectClasses> vote{};
  const auto idx = [](ObjectClass k) { return static_cast<std::size_t>(k); };
  if (looks_static) {
    vote[idx(ObjectClass::kStatic)] = 1.0;
  } else {
    vote[idx(ObjectClass::kStatic)] = -1.0;
    vote[idx(ObjectClass::kVehicle)] =
        (t.rcs_dbsm > cls_cfg_.vehicle_rcs_dbsm ? 1.0 : -0.5) +
        (t.extent_m > cls_cfg_.vehicle_extent_m ? 0.7 : 0.0) +
        (t.velocity_spread_mps < cls_cfg_.pedestrian_spread_mps ? 0.3 : -0.5);
    vote[idx(ObjectClass::kPedestrian)] =
        (t.rcs_dbsm < cls_cfg_.pedestrian_rcs_dbsm ? 1.0 : -1.0) +
        (t.velocity_spread_mps > cls_cfg_.pedestrian_spread_mps ? 1.0 : -0.3) +
        (speed < cls_cfg_.pedestrian_max_speed_mps ? 0.3 : -1.5);
  }

  std::size_t best = 0;
  for (std::size_t k = 0; k < t.evidence.size(); ++k) {
    t.evidence[k] = cls_cfg_.evidence_decay * t.evidence[k] + vote[k];
    if (t.evidence[k] > t.evidence[best]) best = k;
  }
  // Softmax over class evidence as a confidence measure.
  double z = 0.0;
  for (double e : t.evidence) z += std::exp(e - t.evidence[best]);
  t.class_confidence = 1.0 / z;
  t.cls = t.evidence[best] >= cls_cfg_.decision_threshold ? static_cast<ObjectClass>(best)
                                                          : ObjectClass::kUnknown;
}

bool MultiTargetTracker::within_object_extent(const Track& t, const Vec2& p) const {
  // Oriented gate: an extended object spreads along its direction of travel
  // (a vehicle's side), but two objects driving side by side are separated
  // laterally. Slow tracks have no reliable heading and use a round gate.
  const Vec2 d = p - t.position();
  const double speed = t.velocity().norm();
  if (speed < 1.0) return d.norm() < cfg_.absorb_lateral_m;
  const Vec2 along = t.velocity() / speed;
  const double lon = d.dot(along);
  const double lat = d.x() * along.y() - d.y() * along.x();
  return std::abs(lon) < cfg_.absorb_longitudinal_m && std::abs(lat) < cfg_.absorb_lateral_m;
}

void MultiTargetTracker::merge_duplicates() {
  // An extended object occasionally splits into two clusters and spawns a
  // second track; keep the older of two confirmed tracks on the same object.
  std::vector<bool> drop(tracks_.size(), false);
  for (std::size_t i = 0; i < tracks_.size(); ++i) {
    if (!tracks_[i].confirmed() || drop[i]) continue;
    for (std::size_t j = i + 1; j < tracks_.size(); ++j) {
      if (!tracks_[j].confirmed() || drop[j]) continue;
      const Track& older = tracks_[i].age >= tracks_[j].age ? tracks_[i] : tracks_[j];
      const Track& newer = &older == &tracks_[i] ? tracks_[j] : tracks_[i];
      const Vec2 dv = older.velocity() - newer.velocity();
      const Mat2 Sv = older.P.bottomRightCorner<2, 2>() + newer.P.bottomRightCorner<2, 2>();
      const bool same_velocity = dv.norm() < cfg_.merge_velocity_mps ||
                                 dv.dot(Sv.ldlt().solve(dv)) < cfg_.merge_velocity_chi2;
      if (within_object_extent(older, newer.position()) && same_velocity) {
        drop[&older == &tracks_[i] ? j : i] = true;
      }
    }
  }
  std::size_t w = 0;
  for (std::size_t i = 0; i < tracks_.size(); ++i) {
    if (!drop[i]) tracks_[w++] = std::move(tracks_[i]);
  }
  tracks_.resize(w);
}

void MultiTargetTracker::update(double time_s, const std::vector<Cluster>& clusters,
                                const SensorState& sensor) {
  const double dt = initialised_ ? time_s - last_time_s_ : 0.0;
  last_time_s_ = time_s;
  initialised_ = true;
  for (Track& t : tracks_) predict(t, dt);

  std::vector<Measurement> meas;
  meas.reserve(clusters.size());
  for (const Cluster& c : clusters) meas.push_back(make_measurement(c, sensor));

  // Gated cost: Mahalanobis distance + log|S| (negative log-likelihood up to a constant).
  Eigen::MatrixXd cost =
      Eigen::MatrixXd::Constant(static_cast<Eigen::Index>(tracks_.size()),
                                static_cast<Eigen::Index>(meas.size()), kForbidden);
  for (std::size_t i = 0; i < tracks_.size(); ++i) {
    Eigen::Matrix<double, 3, 4> H;
    const Eigen::Vector3d zp = h(tracks_[i], sensor, &H);
    const Eigen::Matrix4d& P = tracks_[i].P;
    for (std::size_t j = 0; j < meas.size(); ++j) {
      const Eigen::Matrix3d S = H * P * H.transpose() + meas[j].R;
      const Eigen::LDLT<Eigen::Matrix3d> ldlt(S);
      const Eigen::Vector3d y = meas[j].z - zp;
      const double d2 = y.dot(ldlt.solve(y));
      if (d2 <= cfg_.gate_chi2) {
        cost(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
            d2 + std::log(ldlt.vectorD().prod());
      }
    }
  }
  const std::vector<int> assignment = solve_assignment(cost, kForbidden);

  std::vector<bool> used(meas.size(), false);
  for (std::size_t i = 0; i < tracks_.size(); ++i) {
    Track& t = tracks_[i];
    ++t.age;
    const int j = assignment[i];
    if (j >= 0) {
      const auto sj = static_cast<std::size_t>(j);
      used[sj] = true;
      correct(t, meas[sj], sensor);
      classify(t, clusters[sj], sensor, time_s);
      ++t.hits;
      t.consecutive_misses = 0;
      t.last_update_s = time_s;
      if (!t.confirmed() && t.hits >= cfg_.confirm_hits) t.status = TrackStatus::kConfirmed;
    } else {
      ++t.consecutive_misses;
    }
  }

  tracks_.erase(std::remove_if(tracks_.begin(), tracks_.end(),
                               [&](const Track& t) {
                                 const int limit = t.confirmed() ? cfg_.max_confirmed_misses
                                                                 : cfg_.max_tentative_misses;
                                 return t.consecutive_misses > limit ||
                                        t.P.topLeftCorner<2, 2>().trace() > 400.0;
                               }),
                tracks_.end());

  for (std::size_t j = 0; j < meas.size(); ++j) {
    if (used[j]) continue;
    // A leftover cluster on the body of an existing confirmed track (another
    // part of the same vehicle) must not start a new track.
    bool absorbed = false;
    for (const Track& t : tracks_) {
      if (!t.confirmed()) continue;
      const double vr_pred = h(t, sensor, nullptr)(2);
      if (std::abs(meas[j].z(2) - vr_pred) < cfg_.absorb_velocity_mps &&
          within_object_extent(t, meas[j].z.head<2>())) {
        absorbed = true;
        break;
      }
    }
    if (absorbed) continue;
    Track t = initiate(meas[j], sensor, time_s);
    classify(t, clusters[j], sensor, time_s);
    tracks_.push_back(std::move(t));
  }
  merge_duplicates();
}

}  // namespace radar
