#pragma once

#include <cstddef>
#include <utility>

#include "radar/types.hpp"

namespace radar::test {

inline constexpr double kDeg = 3.14159265358979323846 / 180.0;

/// (range_bin, doppler_bin) of the strongest cell.
inline std::pair<std::size_t, std::size_t> argmax(const Grid2D<float>& g) {
  std::size_t br = 0, bd = 0;
  for (std::size_t r = 0; r < g.rows(); ++r) {
    for (std::size_t d = 0; d < g.cols(); ++d) {
      if (g(r, d) > g(br, bd)) {
        br = r;
        bd = d;
      }
    }
  }
  return {br, bd};
}

}  // namespace radar::test
