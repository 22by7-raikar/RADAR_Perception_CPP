#include "radar/assignment.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace radar {

std::vector<int> solve_assignment(const Eigen::MatrixXd& cost, double forbidden) {
  const auto rows = static_cast<int>(cost.rows());
  const auto cols = static_cast<int>(cost.cols());
  std::vector<int> result(static_cast<std::size_t>(rows), -1);
  if (rows == 0 || cols == 0) return result;

  // Square-pad; forbidden pairs get a large but finite cost so a perfect
  // matching always exists, and are dropped afterwards.
  const int n = std::max(rows, cols);
  double max_valid = 0.0;
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (cost(i, j) < forbidden) max_valid = std::max(max_valid, std::abs(cost(i, j)));
    }
  }
  const double big = (max_valid + 1.0) * static_cast<double>(n + 1);
  auto c = [&](int i, int j) {
    if (i >= rows || j >= cols) return big;
    const double v = cost(i, j);
    return v >= forbidden ? big : v;
  };

  // Shortest augmenting path formulation with potentials (1-indexed).
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> u(static_cast<std::size_t>(n + 1), 0.0),
      v(static_cast<std::size_t>(n + 1), 0.0);
  std::vector<int> p(static_cast<std::size_t>(n + 1), 0), way(static_cast<std::size_t>(n + 1), 0);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> minv(static_cast<std::size_t>(n + 1), inf);
    std::vector<char> used(static_cast<std::size_t>(n + 1), 0);
    do {
      used[static_cast<std::size_t>(j0)] = 1;
      const int i0 = p[static_cast<std::size_t>(j0)];
      double delta = inf;
      int j1 = 0;
      for (int j = 1; j <= n; ++j) {
        if (used[static_cast<std::size_t>(j)]) continue;
        const double cur =
            c(i0 - 1, j - 1) - u[static_cast<std::size_t>(i0)] - v[static_cast<std::size_t>(j)];
        if (cur < minv[static_cast<std::size_t>(j)]) {
          minv[static_cast<std::size_t>(j)] = cur;
          way[static_cast<std::size_t>(j)] = j0;
        }
        if (minv[static_cast<std::size_t>(j)] < delta) {
          delta = minv[static_cast<std::size_t>(j)];
          j1 = j;
        }
      }
      for (int j = 0; j <= n; ++j) {
        const auto sj = static_cast<std::size_t>(j);
        if (used[sj]) {
          u[static_cast<std::size_t>(p[sj])] += delta;
          v[sj] -= delta;
        } else {
          minv[sj] -= delta;
        }
      }
      j0 = j1;
    } while (p[static_cast<std::size_t>(j0)] != 0);
    do {
      const int j1 = way[static_cast<std::size_t>(j0)];
      p[static_cast<std::size_t>(j0)] = p[static_cast<std::size_t>(j1)];
      j0 = j1;
    } while (j0 != 0);
  }

  for (int j = 1; j <= n; ++j) {
    const int i = p[static_cast<std::size_t>(j)] - 1;
    if (i < rows && j - 1 < cols && cost(i, j - 1) < forbidden) {
      result[static_cast<std::size_t>(i)] = j - 1;
    }
  }
  return result;
}

}  // namespace radar
