#include "radar/clustering.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>

namespace radar {

bool Dbscan::neighbours(const Detection& a, const Detection& b) const {
  const double dx = static_cast<double>(a.x_m - b.x_m);
  const double dy = static_cast<double>(a.y_m - b.y_m);
  const double dv = static_cast<double>(a.velocity_mps - b.velocity_mps);
  const double ep2 = cfg_.eps_position_m * cfg_.eps_position_m;
  const double ev2 = cfg_.eps_velocity_mps * cfg_.eps_velocity_mps;
  return (dx * dx + dy * dy) / ep2 + dv * dv / ev2 <= 1.0;
}

void Dbscan::cluster(const std::vector<Detection>& dets, std::vector<Cluster>& out) {
  out.clear();
  const std::size_t n = dets.size();
  constexpr int kUnvisited = -2, kNoise = -1;
  labels_.assign(n, kUnvisited);

  auto region = [&](std::size_t i, std::vector<std::size_t>& nb) {
    nb.clear();
    for (std::size_t j = 0; j < n; ++j) {
      if (neighbours(dets[i], dets[j])) nb.push_back(j);  // includes i itself
    }
  };

  int next_label = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (labels_[i] != kUnvisited) continue;
    region(i, scratch_);
    if (scratch_.size() < cfg_.min_points) {
      labels_[i] = kNoise;
      continue;
    }
    const int label = next_label++;
    labels_[i] = label;
    queue_.assign(scratch_.begin(), scratch_.end());
    for (std::size_t q = 0; q < queue_.size(); ++q) {
      const std::size_t j = queue_[q];
      if (labels_[j] == kNoise) labels_[j] = label;  // border point
      if (labels_[j] != kUnvisited) continue;
      labels_[j] = label;
      region(j, scratch_);
      if (scratch_.size() >= cfg_.min_points) {
        queue_.insert(queue_.end(), scratch_.begin(), scratch_.end());
      }
    }
  }

  out.resize(static_cast<std::size_t>(next_label));
  for (std::size_t i = 0; i < n; ++i) {
    if (labels_[i] >= 0) out[static_cast<std::size_t>(labels_[i])].members.push_back(i);
  }

  for (Cluster& c : out) {
    double w_sum = 0.0, v_min = std::numeric_limits<double>::infinity(), v_max = -v_min;
    Vec2 mean = Vec2::Zero();
    double v_mean = 0.0, snr_max = -std::numeric_limits<double>::infinity();
    for (std::size_t i : c.members) {
      const Detection& d = dets[i];
      const double w = std::pow(10.0, static_cast<double>(d.rcs_dbsm) / 10.0);
      mean += w * Vec2(d.x_m, d.y_m);
      v_mean += w * static_cast<double>(d.velocity_mps);
      w_sum += w;
      v_min = std::min(v_min, static_cast<double>(d.velocity_mps));
      v_max = std::max(v_max, static_cast<double>(d.velocity_mps));
      snr_max = std::max(snr_max, static_cast<double>(d.snr_db));
    }
    c.centroid = mean / w_sum;
    c.radial_velocity = v_mean / w_sum;
    c.velocity_spread = v_max - v_min;
    c.rcs_dbsm = 10.0 * std::log10(w_sum);
    c.max_snr_db = snr_max;

    // Extent: full spread along the principal axes of the (unweighted) points.
    if (c.members.size() >= 2) {
      Vec2 m = Vec2::Zero();
      for (std::size_t i : c.members) m += Vec2(dets[i].x_m, dets[i].y_m);
      m /= static_cast<double>(c.members.size());
      Mat2 cov = Mat2::Zero();
      for (std::size_t i : c.members) {
        const Vec2 e = Vec2(dets[i].x_m, dets[i].y_m) - m;
        cov += e * e.transpose();
      }
      Eigen::SelfAdjointEigenSolver<Mat2> es(cov);
      const Vec2 major = es.eigenvectors().col(1), minor = es.eigenvectors().col(0);
      double lo1 = 0, hi1 = 0, lo2 = 0, hi2 = 0;
      for (std::size_t i : c.members) {
        const Vec2 e = Vec2(dets[i].x_m, dets[i].y_m) - m;
        lo1 = std::min(lo1, e.dot(major));
        hi1 = std::max(hi1, e.dot(major));
        lo2 = std::min(lo2, e.dot(minor));
        hi2 = std::max(hi2, e.dot(minor));
      }
      c.extent_major = hi1 - lo1;
      c.extent_minor = hi2 - lo2;
    }
  }
}

}  // namespace radar
