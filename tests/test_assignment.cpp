#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <random>

#include "radar/assignment.hpp"

namespace radar {

namespace {

double brute_force_cost(const Eigen::MatrixXd& c) {
  // Square, all valid: try every permutation.
  std::vector<int> perm(static_cast<std::size_t>(c.cols()));
  std::iota(perm.begin(), perm.end(), 0);
  double best = 1e18;
  do {
    double s = 0.0;
    for (Eigen::Index i = 0; i < c.rows(); ++i) s += c(i, perm[static_cast<std::size_t>(i)]);
    best = std::min(best, s);
  } while (std::next_permutation(perm.begin(), perm.end()));
  return best;
}

}  // namespace

TEST(Assignment, KnownOptimum) {
  Eigen::MatrixXd c(3, 3);
  c << 4, 1, 3, 2, 0, 5, 3, 2, 2;
  const auto a = solve_assignment(c);
  EXPECT_EQ(a, (std::vector<int>{1, 0, 2}));  // 1 + 2 + 2 = 5
}

TEST(Assignment, MatchesBruteForceOnRandomMatrices) {
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> u(0.0, 10.0);
  for (int trial = 0; trial < 50; ++trial) {
    Eigen::MatrixXd c(6, 6);
    for (Eigen::Index i = 0; i < c.size(); ++i) c.data()[i] = u(rng);
    const auto a = solve_assignment(c);
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += c(static_cast<Eigen::Index>(i), a[i]);
    EXPECT_NEAR(s, brute_force_cost(c), 1e-9);
  }
}

TEST(Assignment, RectangularAndForbidden) {
  Eigen::MatrixXd c(2, 3);
  c << 1e9, 5, 1, 2, 1e9, 1e9;
  const auto a = solve_assignment(c, 1e9);
  EXPECT_EQ(a[0], 2);
  EXPECT_EQ(a[1], 0);

  Eigen::MatrixXd all_forbidden = Eigen::MatrixXd::Constant(3, 2, 1e9);
  for (int j : solve_assignment(all_forbidden, 1e9)) EXPECT_EQ(j, -1);
  EXPECT_TRUE(solve_assignment(Eigen::MatrixXd(0, 4)).empty());
}

}  // namespace radar
