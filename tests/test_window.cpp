#include <gtest/gtest.h>

#include "radar/window.hpp"

namespace radar {

TEST(Window, RectangularIsAllOnes) {
  const auto w = make_window(WindowType::kRectangular, 16);
  for (float v : w) EXPECT_FLOAT_EQ(v, 1.0f);
  EXPECT_DOUBLE_EQ(coherent_gain(w), 1.0);
}

TEST(Window, HannIsPeriodicAndHasHalfGain) {
  const std::size_t n = 64;
  const auto w = make_window(WindowType::kHann, n);
  ASSERT_EQ(w.size(), n);
  EXPECT_NEAR(w[0], 0.0f, 1e-7);
  EXPECT_NEAR(w[n / 2], 1.0f, 1e-7);
  for (std::size_t i = 1; i < n; ++i) EXPECT_NEAR(w[i], w[n - i], 1e-6);  // periodic symmetry
  EXPECT_NEAR(coherent_gain(w), 0.5, 1e-6);
}

TEST(Window, KnownCoherentGains) {
  EXPECT_NEAR(coherent_gain(make_window(WindowType::kHamming, 128)), 0.54, 1e-6);
  EXPECT_NEAR(coherent_gain(make_window(WindowType::kBlackman, 128)), 0.42, 1e-6);
}

}  // namespace radar
