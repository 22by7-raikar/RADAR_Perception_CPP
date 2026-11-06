#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>

#include "radar/impairments.hpp"
#include "radar/mitigation.hpp"
#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"
#include "test_helpers.hpp"

namespace radar {

namespace {

using test::kDeg;

Grid2D<float> gamma_noise(std::size_t rows, std::size_t cols, std::size_t looks,
                          std::mt19937& rng) {
  std::gamma_distribution<float> dist(static_cast<float>(looks), 1.0f);
  Grid2D<float> g(rows, cols);
  for (std::size_t i = 0; i < g.size(); ++i) g.data()[i] = dist(rng);
  return g;
}

bool has_cell(const std::vector<CfarCell>& cells, std::size_t r, std::size_t d) {
  return std::any_of(cells.begin(), cells.end(),
                     [&](const CfarCell& c) { return c.range_bin == r && c.doppler_bin == d; });
}

Interferer default_interferer() {
  Interferer in;
  in.time_offset_s = 3e-6;
  in.azimuth_rad = 10.0 * kDeg;
  return in;
}

}  // namespace

// --- OS-CFAR -----------------------------------------------------------------

TEST(OsCfar, SingleLookMatchesClosedForm) {
  // Exponential noise: Pfa = prod_{i<k} (N - i) / (N - i + T).
  const std::size_t n = 24, k = 18;
  const double t = 3.0;
  double expected = 1.0;
  for (std::size_t i = 0; i < k; ++i) {
    expected *= static_cast<double>(n - i) / (static_cast<double>(n - i) + t);
  }
  EXPECT_NEAR(OsCfar2D::false_alarm_probability(t, n, k, 1) / expected, 1.0, 1e-3);
}

TEST(OsCfar, ThresholdInvertsFalseAlarmProbability) {
  const double t = OsCfar2D::threshold_multiplier(120, 90, 8, 1e-5);
  EXPECT_NEAR(OsCfar2D::false_alarm_probability(t, 120, 90, 8) / 1e-5, 1.0, 1e-3);
}

TEST(OsCfar, EmpiricalFalseAlarmRateMatchesDesign) {
  std::mt19937 rng(77);
  CfarConfig c;
  c.pfa = 1e-3;
  c.num_looks = 8;
  OsCfar2D cfar(c);
  std::size_t alarms = 0, tested = 0;
  std::vector<CfarCell> cells;
  for (int trial = 0; trial < 12; ++trial) {
    const auto p = gamma_noise(256, 128, c.num_looks, rng);
    cfar.detect(p, cells);
    alarms += cells.size();
    tested += p.size();
  }
  const double rate = static_cast<double>(alarms) / static_cast<double>(tested);
  EXPECT_GT(rate, 0.8 * c.pfa);
  EXPECT_LT(rate, 1.2 * c.pfa);
}

TEST(OsCfar, WeakTargetBesideStrongOneIsNotMasked) {
  // A 45 dB reflector 4 range bins (inside the training window) from a
  // 16 dB one: the cell-averaging estimate is dominated by the strong cell.
  std::mt19937 rng(3);
  CfarConfig c;
  c.num_looks = 8;
  auto p = gamma_noise(128, 64, c.num_looks, rng);
  p(60, 30) = 8.0f * std::pow(10.0f, 4.5f);
  p(64, 30) = 8.0f * std::pow(10.0f, 1.6f);
  std::vector<CfarCell> ca_cells, os_cells;
  CaCfar2D(c).detect(p, ca_cells);
  OsCfar2D(c).detect(p, os_cells);
  EXPECT_TRUE(has_cell(ca_cells, 60, 30));
  EXPECT_FALSE(has_cell(ca_cells, 64, 30));
  EXPECT_TRUE(has_cell(os_cells, 60, 30));
  EXPECT_TRUE(has_cell(os_cells, 64, 30));
}

// --- Interference ------------------------------------------------------------

TEST(Interference, BurstsAreShortAndRaiseTheNoiseFloor) {
  const RadarConfig radar;
  SimulatorConfig sc;
  Simulator sim(radar, sc);
  DataCube clean = sim.simulate_frame({});
  DataCube dirty = clean;
  add_interference(radar, default_interferer(), dirty);

  std::size_t touched = 0;
  for (std::size_t i = 0; i < clean.size(); ++i) touched += clean.data()[i] != dirty.data()[i];
  const double fraction = static_cast<double>(touched) / static_cast<double>(clean.size());
  EXPECT_GT(fraction, 0.005);
  EXPECT_LT(fraction, 0.10);  // short bursts, not a continuous jammer

  RadarPipeline pipe(radar);
  FrameResult a, b;
  pipe.process(clean, a);
  pipe.process(dirty, b);
  EXPECT_GT(10.0 * std::log10(b.noise_floor / a.noise_floor), 10.0);
}

TEST(Interference, BlankingRestoresFloorAndWeakTarget) {
  const RadarConfig radar;
  SimulatorConfig sc;
  sc.seed = 4;
  Simulator sim(radar, sc);
  Target weak;
  weak.range_m = 150.0;
  weak.velocity_mps = -6.0;
  weak.rcs_dbsm = 0.0;  // ~ 10 dB above the floor after integration
  DataCube clean = sim.simulate_frame({weak});
  DataCube dirty = clean;
  add_interference(radar, default_interferer(), dirty);

  PipelineConfig on;
  on.interference.enabled = true;
  RadarPipeline plain(radar), mitigated(radar, on);
  FrameResult ref, hit, fixed;
  plain.process(clean, ref);
  plain.process(dirty, hit);
  mitigated.process(dirty, fixed);

  auto found = [&](const FrameResult& r) {
    return std::any_of(r.detections.begin(), r.detections.end(),
                       [](const Detection& d) { return std::abs(d.range_m - 150.0) < 2.0; });
  };
  ASSERT_TRUE(found(ref));
  EXPECT_FALSE(found(hit));
  EXPECT_TRUE(found(fixed));
  EXPECT_GT(fixed.interference_samples_blanked, 0u);
  EXPECT_LT(10.0 * std::log10(fixed.noise_floor / ref.noise_floor), 1.0);
}

TEST(Interference, BlankingLeavesCleanMultiTargetFramesAlone) {
  const RadarConfig radar;
  Simulator sim(radar);
  std::vector<Target> targets;
  for (int i = 0; i < 8; ++i) {
    targets.push_back({12.0 + 23.0 * i, -20.0 + 5.5 * i, (-40.0 + 11.0 * i) * kDeg, 15.0});
  }
  PipelineConfig on;
  on.interference.enabled = true;
  RadarPipeline mitigated(radar, on), plain(radar);
  FrameResult r, ref;
  const DataCube frame = sim.simulate_frame(targets);
  mitigated.process(frame, r);
  plain.process(frame, ref);
  // Eight strong echoes must not look like bursts. An occasional single
  // noise trigger (one sample plus guards, then refilled) is harmless.
  EXPECT_LE(r.interference_samples_blanked, 10u);
  EXPECT_EQ(r.detections.size(), ref.detections.size());
  EXPECT_NEAR(r.noise_floor / ref.noise_floor, 1.0, 0.01);
}

TEST(Interference, MimoBlanksEachTransmitSlotSeparately) {
  const RadarConfig radar = RadarConfig::tdm_mimo_3tx4rx();
  SimulatorConfig sc;
  Simulator sim(radar, sc);
  DataCube dirty = sim.simulate_frame({});
  add_interference(radar, default_interferer(), dirty);
  std::vector<float> scratch;
  InterferenceMitigationConfig cfg;
  cfg.enabled = true;
  cfg.reconstruct_iterations = 0;  // leave the blanked slots at zero to inspect them
  EXPECT_GT(mitigate_interference(radar, cfg, dirty, scratch), 0u);
  // Every blanked slot is zero across the four receivers of its transmitter.
  for (std::size_t g = 0; g < radar.num_tx; ++g) {
    for (std::size_t m = 0; m < radar.chirps_per_frame; ++m) {
      for (std::size_t n = 0; n < radar.samples_per_chirp; ++n) {
        const bool zero0 = dirty(g * 4, m, n) == cfloat{};
        for (std::size_t r = 1; r < 4; ++r) ASSERT_EQ(dirty(g * 4 + r, m, n) == cfloat{}, zero0);
      }
    }
  }
}

// --- Multipath -----------------------------------------------------------------

TEST(Multipath, MirrorGeometry) {
  PlanarReflector wall;
  wall.y_m = -4.0;
  wall.loss_db = 3.0;
  MovingPoint car{{40.0, 0.0}, {-10.0, 0.0}, 10.0};
  const auto ghosts = multipath_targets(car, wall);
  ASSERT_EQ(ghosts.size(), 3u);
  const double r_img = std::hypot(40.0, 8.0);
  EXPECT_NEAR(ghosts[0].range_m, r_img, 1e-9);  // double bounce at the image
  EXPECT_NEAR(ghosts[0].azimuth_rad, std::atan2(-8.0, 40.0), 1e-9);
  EXPECT_NEAR(ghosts[0].rcs_dbsm, 4.0, 1e-9);
  EXPECT_NEAR(ghosts[1].range_m, 0.5 * (40.0 + r_img), 1e-9);  // single bounce
  EXPECT_NEAR(ghosts[1].azimuth_rad, 0.0, 1e-9);
  EXPECT_NEAR(ghosts[2].azimuth_rad, std::atan2(-8.0, 40.0), 1e-9);

  MovingPoint behind{{40.0, -6.0}, {0.0, 0.0}, 10.0};
  EXPECT_TRUE(multipath_targets(behind, wall).empty());
}

namespace {

Detection make_det(double x, double y, double v) {
  Detection d;
  d.x_m = static_cast<float>(x);
  d.y_m = static_cast<float>(y);
  d.range_m = static_cast<float>(std::hypot(x, y));
  d.azimuth_rad = static_cast<float>(std::atan2(y, x));
  d.velocity_mps = static_cast<float>(v);
  return d;
}

}  // namespace

TEST(Multipath, FitsGuardrailAndFlagsMirroredGhost) {
  const Vec2 ego(20.0, 0.0);
  std::vector<Detection> dets;
  std::mt19937 rng(2);
  std::normal_distribution<double> jitter(0.0, 0.1);
  for (double x = 8.0; x < 100.0; x += 4.0) {  // guardrail posts at y = -4
    const double y = -4.0 + jitter(rng);
    dets.push_back(make_det(x, y, -20.0 * x / std::hypot(x, y)));
  }
  dets.push_back(make_det(40.0, 4.5, -5.0));    // tree / pole on the other side
  dets.push_back(make_det(50.0, 0.5, -3.0));    // real car ahead
  dets.push_back(make_det(50.0, -8.5, -3.05));  // its mirror image behind the rail
  dets.push_back(make_det(70.0, -9.0, 4.0));    // real vehicle on a parallel road

  GhostFilterConfig cfg;
  const ReflectorLine right = fit_reflector(dets, ego, -1, cfg);
  ASSERT_TRUE(right.valid);
  EXPECT_NEAR(right.slope, 0.0, 0.01);
  EXPECT_NEAR(right.offset, -4.0, 0.15);
  EXPECT_FALSE(fit_reflector(dets, ego, +1, cfg).valid);  // one point is not a rail

  std::vector<char> ghost;
  flag_multipath_ghosts(dets, {right}, ego, cfg, ghost);
  const std::size_t n = dets.size();
  EXPECT_FALSE(ghost[n - 3]);  // the car
  EXPECT_TRUE(ghost[n - 2]);   // its mirror image
  EXPECT_FALSE(ghost[n - 1]);  // behind the rail but has no mirror twin
  EXPECT_EQ(std::count(ghost.begin(), ghost.end(), 1), 1);
}

}  // namespace radar
