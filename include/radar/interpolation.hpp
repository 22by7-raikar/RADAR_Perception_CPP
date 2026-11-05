#pragma once

#include <algorithm>
#include <cmath>

namespace radar {

/// Sub-bin peak offset in [-0.5, 0.5] from a parabola through three samples.
inline double parabolic_offset(double left, double centre, double right) {
  const double denom = left - 2.0 * centre + right;
  if (denom >= 0.0) return 0.0;  // not a strict local maximum
  return std::clamp(0.5 * (left - right) / denom, -0.5, 0.5);
}

/// Parabolic interpolation on log power, which is exact for Gaussian-shaped
/// peaks and a good fit for Hann-windowed main lobes.
inline double log_parabolic_offset(double left, double centre, double right) {
  constexpr double kFloor = 1e-30;
  return parabolic_offset(std::log(std::max(left, kFloor)), std::log(std::max(centre, kFloor)),
                          std::log(std::max(right, kFloor)));
}

}  // namespace radar
