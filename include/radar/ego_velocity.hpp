#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "radar/geometry.hpp"
#include "radar/types.hpp"

namespace radar {

struct EgoVelocityConfig {
  double inlier_threshold_mps = 0.25;  ///< |measured - predicted| radial velocity
  int iterations = 80;
  std::size_t min_inliers = 6;
  double min_inlier_fraction = 0.3;
  double min_angle_separation_rad = 0.05;  ///< reject near-degenerate 2-point samples
  double residual_floor_mps = 0.08;        ///< lower bound on the noise std in the covariance
  /// Variance inflation: point residuals are not independent (shared angle
  /// bias, cell sharing), so the textbook LS covariance is optimistic. 9x
  /// (3x in std) matches the empirical error in the drive scenario.
  double covariance_inflation = 9.0;
  std::uint32_t seed = 1;
};

/// Sensor velocity over ground, in the sensor frame, from one radar scan.
struct EgoVelocityEstimate {
  bool valid = false;
  Vec2 velocity = Vec2::Zero();
  Mat2 covariance = Mat2::Identity();
  std::size_t inliers = 0;
  double residual_std_mps = 0.0;
  std::vector<char> inlier_mask;  ///< per detection: 1 = consistent with a static world
};

/// Instantaneous ego-velocity from Doppler.
///
/// Every stationary reflector at azimuth theta obeys
///   v_r = -(cos(theta), sin(theta)) . v_sensor,
/// a linear model in the two sensor-velocity components. RANSAC on 2-point
/// minimal samples rejects moving objects, then least squares on the inliers
/// gives the estimate and its covariance.
class EgoVelocityEstimator {
 public:
  explicit EgoVelocityEstimator(const EgoVelocityConfig& cfg = {});

  EgoVelocityEstimate estimate(const std::vector<Detection>& detections);

 private:
  bool least_squares(const std::vector<Detection>& d, const std::vector<char>& mask, Vec2& v,
                     Mat2* ata_inv) const;
  std::size_t count_inliers(const std::vector<Detection>& d, const Vec2& v,
                            std::vector<char>& mask) const;

  EgoVelocityConfig cfg_;
  std::mt19937 rng_;
};

}  // namespace radar
