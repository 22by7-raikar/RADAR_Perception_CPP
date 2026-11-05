#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>

#include "radar/io.hpp"

namespace radar {

TEST(Io, LoadsScenarioWithCommentsAndHeader) {
  const auto path = std::filesystem::temp_directory_path() / "rp_scenario_test.csv";
  {
    std::ofstream f(path);
    f << "# comment\nrange_m,velocity_mps,azimuth_deg,rcs_dbsm\n\n10,-2,90,5\n20,3,0\n";
  }
  const auto t = load_targets_csv(path.string());
  ASSERT_EQ(t.size(), 2u);
  EXPECT_DOUBLE_EQ(t[0].range_m, 10.0);
  EXPECT_DOUBLE_EQ(t[0].velocity_mps, -2.0);
  EXPECT_NEAR(t[0].azimuth_rad, kPi / 2.0, 1e-12);
  EXPECT_DOUBLE_EQ(t[0].rcs_dbsm, 5.0);
  EXPECT_DOUBLE_EQ(t[1].rcs_dbsm, Target{}.rcs_dbsm);
  std::filesystem::remove(path);
}

TEST(Io, ReportsMalformedLines) {
  const auto path = std::filesystem::temp_directory_path() / "rp_scenario_bad.csv";
  {
    std::ofstream f(path);
    f << "10,abc,0\n";
  }
  EXPECT_THROW(load_targets_csv(path.string()), std::runtime_error);
  std::filesystem::remove(path);
}

TEST(Io, MissingFileThrows) {
  EXPECT_THROW(load_targets_csv("/nonexistent/scenario.csv"), std::runtime_error);
}

}  // namespace radar
