#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <set>
#include <utility>

#include "radar/cfar.hpp"

namespace radar {

namespace {

Grid2D<float> gamma_noise(std::size_t rows, std::size_t cols, std::size_t looks,
                          std::mt19937& rng) {
  std::gamma_distribution<float> dist(static_cast<float>(looks), 1.0f);
  Grid2D<float> g(rows, cols);
  for (std::size_t i = 0; i < g.size(); ++i) g.data()[i] = dist(rng);
  return g;
}

// Reference implementation: explicit loops over the window, same edge rules.
std::set<std::pair<std::size_t, std::size_t>> brute_force(const Grid2D<float>& p,
                                                          const CfarConfig& c) {
  std::set<std::pair<std::size_t, std::size_t>> hits;
  const auto nr = static_cast<long>(p.rows());
  const auto nd = static_cast<long>(p.cols());
  const long wr = static_cast<long>(c.train_range + c.guard_range);
  const long wd = static_cast<long>(c.train_doppler + c.guard_doppler);
  const long gr = static_cast<long>(c.guard_range);
  const long gd = static_cast<long>(c.guard_doppler);
  for (long r = 0; r < nr; ++r) {
    for (long d = 0; d < nd; ++d) {
      double sum = 0.0;
      std::size_t n = 0;
      for (long i = r - wr; i <= r + wr; ++i) {
        if (i < 0 || i >= nr) continue;
        for (long j = d - wd; j <= d + wd; ++j) {
          if (std::abs(i - r) <= gr && std::abs(j - d) <= gd) continue;
          sum += p(static_cast<std::size_t>(i), static_cast<std::size_t>((j + nd) % nd));
          ++n;
        }
      }
      const double t = CaCfar2D::threshold_multiplier(n, c.num_looks, c.pfa);
      if (p(static_cast<std::size_t>(r), static_cast<std::size_t>(d)) > t * sum) {
        hits.insert({static_cast<std::size_t>(r), static_cast<std::size_t>(d)});
      }
    }
  }
  return hits;
}

}  // namespace

TEST(CfarThreshold, SingleLookMatchesClosedForm) {
  for (std::size_t n : {8u, 24u, 100u}) {
    for (double pfa : {1e-3, 1e-6}) {
      const double expected = std::pow(pfa, -1.0 / static_cast<double>(n)) - 1.0;
      EXPECT_NEAR(CaCfar2D::threshold_multiplier(n, 1, pfa), expected, 1e-9 * expected);
    }
  }
}

TEST(CfarThreshold, InvertsFalseAlarmProbability) {
  for (std::size_t looks : {1u, 4u, 8u}) {
    const double t = CaCfar2D::threshold_multiplier(48, looks, 1e-5);
    EXPECT_NEAR(CaCfar2D::false_alarm_probability(t, 48, looks), 1e-5, 1e-9);
  }
}

TEST(Cfar, SummedAreaTableMatchesBruteForce) {
  std::mt19937 rng(7);
  CfarConfig c;
  c.train_range = 4;
  c.train_doppler = 3;
  c.guard_range = 1;
  c.guard_doppler = 2;
  c.pfa = 1e-2;
  c.num_looks = 2;
  auto p = gamma_noise(48, 32, c.num_looks, rng);
  p(0, 0) = 60.0f;    // corner: clipped range + wrapped Doppler
  p(20, 31) = 60.0f;  // Doppler edge
  p(47, 10) = 60.0f;  // range edge

  CaCfar2D cfar(c);
  std::vector<CfarCell> cells;
  cfar.detect(p, cells);
  std::set<std::pair<std::size_t, std::size_t>> fast;
  for (const auto& cell : cells) fast.insert({cell.range_bin, cell.doppler_bin});

  EXPECT_EQ(fast, brute_force(p, c));
  EXPECT_TRUE(fast.count({0, 0}));
  EXPECT_TRUE(fast.count({20, 31}));
  EXPECT_TRUE(fast.count({47, 10}));
}

TEST(Cfar, EmpiricalFalseAlarmRateMatchesDesign) {
  std::mt19937 rng(1234);
  CfarConfig c;
  c.pfa = 1e-3;
  c.num_looks = 8;
  CaCfar2D cfar(c);

  std::size_t alarms = 0, cells_tested = 0;
  std::vector<CfarCell> cells;
  for (int trial = 0; trial < 20; ++trial) {
    const auto p = gamma_noise(256, 128, c.num_looks, rng);
    cfar.detect(p, cells);
    alarms += cells.size();
    cells_tested += p.size();
  }
  const double rate = static_cast<double>(alarms) / static_cast<double>(cells_tested);
  EXPECT_GT(rate, 0.8 * c.pfa);
  EXPECT_LT(rate, 1.2 * c.pfa);
}

TEST(Cfar, NoiseEstimateIsTrainingMean) {
  Grid2D<float> p(32, 32, 1.0f);
  p(16, 16) = 1000.0f;
  CfarConfig c;
  c.num_looks = 1;
  CaCfar2D cfar(c);
  std::vector<CfarCell> cells;
  cfar.detect(p, cells);
  ASSERT_EQ(cells.size(), 1u);
  EXPECT_EQ(cells[0].range_bin, 16u);
  EXPECT_EQ(cells[0].doppler_bin, 16u);
  EXPECT_NEAR(cells[0].noise, 1.0f, 1e-5);
}

TEST(Cfar, RejectsOversizedWindow) {
  Grid2D<float> p(8, 8, 1.0f);
  CaCfar2D cfar(CfarConfig{});
  std::vector<CfarCell> cells;
  EXPECT_THROW(cfar.detect(p, cells), std::invalid_argument);
}

TEST(PeakGrouping, CollapsesBlobToSinglePeak) {
  Grid2D<float> p(16, 16, 1.0f);
  std::vector<CfarCell> cells;
  // Main-lobe-like blob: power falls off with distance from the peak at (7, 8).
  for (std::size_t r = 5; r <= 9; ++r) {
    for (std::size_t d = 5; d <= 10; ++d) {
      const float dist =
          std::abs(static_cast<float>(r) - 7.0f) + std::abs(static_cast<float>(d) - 8.0f);
      p(r, d) = 80.0f - 10.0f * dist;
      cells.push_back({r, d, p(r, d), 1.0f});
    }
  }
  std::vector<CfarCell> peaks;
  suppress_non_maxima(p, cells, peaks);
  ASSERT_EQ(peaks.size(), 1u);
  EXPECT_EQ(peaks[0].range_bin, 7u);
  EXPECT_EQ(peaks[0].doppler_bin, 8u);
}

TEST(PeakGrouping, FlatTopYieldsOnePeak) {
  Grid2D<float> p(8, 8, 1.0f);
  p(3, 3) = 10.0f;
  p(3, 4) = 10.0f;
  std::vector<CfarCell> cells{{3, 3, 10.0f, 1.0f}, {3, 4, 10.0f, 1.0f}};
  std::vector<CfarCell> peaks;
  suppress_non_maxima(p, cells, peaks);
  EXPECT_EQ(peaks.size(), 1u);
}

}  // namespace radar
