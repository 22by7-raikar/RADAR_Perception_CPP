#pragma once

#include <Eigen/Core>
#include <vector>

namespace radar {

/// Minimum-cost assignment (Hungarian / Kuhn-Munkres, O(n^3)) for a
/// rectangular cost matrix. Entries >= `forbidden` are never assigned.
///
/// Returns, for each row, the assigned column or -1.
std::vector<int> solve_assignment(const Eigen::MatrixXd& cost, double forbidden = 1e9);

}  // namespace radar
