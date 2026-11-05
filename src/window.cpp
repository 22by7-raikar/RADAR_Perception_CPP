#include "radar/window.hpp"

#include <cmath>
#include <numeric>

#include "radar/config.hpp"

namespace radar {

std::vector<float> make_window(WindowType type, std::size_t n) {
  std::vector<float> w(n, 1.0f);
  if (n <= 1 || type == WindowType::kRectangular) return w;

  const double step = 2.0 * kPi / static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = step * static_cast<double>(i);
    double v = 1.0;
    switch (type) {
      case WindowType::kHann:
        v = 0.5 - 0.5 * std::cos(x);
        break;
      case WindowType::kHamming:
        v = 0.54 - 0.46 * std::cos(x);
        break;
      case WindowType::kBlackman:
        v = 0.42 - 0.5 * std::cos(x) + 0.08 * std::cos(2.0 * x);
        break;
      case WindowType::kRectangular:
        break;
    }
    w[i] = static_cast<float>(v);
  }
  return w;
}

double coherent_gain(const std::vector<float>& window) {
  if (window.empty()) return 0.0;
  return std::accumulate(window.begin(), window.end(), 0.0) / static_cast<double>(window.size());
}

}  // namespace radar
