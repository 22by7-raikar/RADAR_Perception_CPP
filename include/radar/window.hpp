#pragma once

#include <cstddef>
#include <vector>

namespace radar {

enum class WindowType { kRectangular, kHann, kHamming, kBlackman };

/// Periodic (DFT-even) taper of length n.
std::vector<float> make_window(WindowType type, std::size_t n);

/// Mean of the window coefficients (amplitude gain for a bin-centred tone).
double coherent_gain(const std::vector<float>& window);

}  // namespace radar
