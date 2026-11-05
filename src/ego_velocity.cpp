#include "radar/ego_velocity.hpp"

#include <Eigen/LU>
#include <cmath>

namespace radar {

namespace {

Vec2 row(const Detection& d) {
  const double az = static_cast<double>(d.azimuth_rad);
  return {-std::cos(az), -std::sin(az)};
}

}  // namespace

EgoVelocityEstimator::EgoVelocityEstimator(const EgoVelocityConfig& cfg)
    : cfg_(cfg), rng_(cfg.seed) {}

std::size_t EgoVelocityEstimator::count_inliers(const std::vector<Detection>& d, const Vec2& v,
                                                std::vector<char>& mask) const {
  std::size_t n = 0;
  mask.assign(d.size(), 0);
  for (std::size_t i = 0; i < d.size(); ++i) {
    const double r = static_cast<double>(d[i].velocity_mps) - row(d[i]).dot(v);
    if (std::abs(r) < cfg_.inlier_threshold_mps) {
      mask[i] = 1;
      ++n;
    }
  }
  return n;
}

bool EgoVelocityEstimator::least_squares(const std::vector<Detection>& d,
                                         const std::vector<char>& mask, Vec2& v,
                                         Mat2* ata_inv) const {
  Mat2 ata = Mat2::Zero();
  Vec2 atb = Vec2::Zero();
  for (std::size_t i = 0; i < d.size(); ++i) {
    if (!mask[i]) continue;
    const Vec2 a = row(d[i]);
    ata += a * a.transpose();
    atb += a * static_cast<double>(d[i].velocity_mps);
  }
  if (std::abs(ata.determinant()) < 1e-9) return false;
  const Mat2 inv = ata.inverse();
  v = inv * atb;
  if (ata_inv) *ata_inv = inv;
  return true;
}

EgoVelocityEstimate EgoVelocityEstimator::estimate(const std::vector<Detection>& d) {
  EgoVelocityEstimate est;
  est.inlier_mask.assign(d.size(), 0);
  if (d.size() < cfg_.min_inliers) return est;

  std::uniform_int_distribution<std::size_t> pick(0, d.size() - 1);
  std::vector<char> mask;
  std::size_t best_count = 0;
  Vec2 best_v = Vec2::Zero();
  for (int it = 0; it < cfg_.iterations; ++it) {
    const std::size_t i = pick(rng_), j = pick(rng_);
    if (i == j) continue;
    if (std::abs(static_cast<double>(d[i].azimuth_rad - d[j].azimuth_rad)) <
        cfg_.min_angle_separation_rad) {
      continue;
    }
    Mat2 A;
    A.row(0) = row(d[i]).transpose();
    A.row(1) = row(d[j]).transpose();
    const Vec2 b(d[i].velocity_mps, d[j].velocity_mps);
    if (std::abs(A.determinant()) < 1e-6) continue;
    const Vec2 v = A.inverse() * b;
    const std::size_t n = count_inliers(d, v, mask);
    if (n > best_count) {
      best_count = n;
      best_v = v;
    }
  }
  if (best_count < cfg_.min_inliers ||
      static_cast<double>(best_count) < cfg_.min_inlier_fraction * static_cast<double>(d.size())) {
    return est;
  }

  // Refine: least squares on the consensus set, then re-select inliers once.
  Vec2 v = best_v;
  Mat2 ata_inv;
  count_inliers(d, v, mask);
  if (!least_squares(d, mask, v, &ata_inv)) return est;
  const std::size_t n = count_inliers(d, v, mask);
  if (n < cfg_.min_inliers || !least_squares(d, mask, v, &ata_inv)) return est;

  double ss = 0.0;
  for (std::size_t i = 0; i < d.size(); ++i) {
    if (!mask[i]) continue;
    const double r = static_cast<double>(d[i].velocity_mps) - row(d[i]).dot(v);
    ss += r * r;
  }
  const double dof = std::max(1.0, static_cast<double>(n) - 2.0);
  const double sigma = std::max(cfg_.residual_floor_mps, std::sqrt(ss / dof));

  est.valid = true;
  est.velocity = v;
  est.covariance = cfg_.covariance_inflation * sigma * sigma * ata_inv;
  est.inliers = n;
  est.residual_std_mps = sigma;
  est.inlier_mask = mask;
  return est;
}

}  // namespace radar
