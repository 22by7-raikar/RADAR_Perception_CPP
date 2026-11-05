#include <gtest/gtest.h>

#include <cmath>

#include "radar/evaluation.hpp"
#include "radar/io.hpp"
#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"

namespace radar {

namespace {

FrameResult run(const std::vector<Target>& targets, std::uint32_t seed = 42) {
  RadarConfig cfg;
  SimulatorConfig sc;
  sc.seed = seed;
  Simulator sim(cfg, sc);
  RadarPipeline pipe(cfg);
  FrameResult res;
  pipe.process(sim.simulate_frame(targets), res);
  return res;
}

void expect_all_found(const std::vector<Target>& truth, const FrameResult& res,
                      std::size_t max_false_alarms) {
  const auto ev = evaluate(truth, res.detections);
  EXPECT_EQ(ev.false_negatives, 0u);
  EXPECT_LE(ev.false_positives, max_false_alarms);
  for (std::size_t i = 0; i < truth.size(); ++i) {
    ASSERT_GE(ev.truth_to_detection[i], 0) << "target " << i << " missed";
    const Detection& d = res.detections[static_cast<std::size_t>(ev.truth_to_detection[i])];
    EXPECT_NEAR(d.range_m, truth[i].range_m, 0.5) << "target " << i;
    EXPECT_NEAR(d.velocity_mps, truth[i].velocity_mps, 0.5) << "target " << i;
    EXPECT_NEAR(d.azimuth_rad, truth[i].azimuth_rad, 3.0 * kPi / 180.0) << "target " << i;
    EXPECT_NEAR(d.x_m, truth[i].range_m * std::cos(truth[i].azimuth_rad), 3.0);
    EXPECT_NEAR(d.y_m, truth[i].range_m * std::sin(truth[i].azimuth_rad), 3.0);
  }
}

}  // namespace

TEST(Pipeline, DetectsAllTargetsInMultiTargetScene) {
  const auto truth = load_targets_csv(std::string(RP_SCENARIO_DIR) + "/multi_target.csv");
  ASSERT_EQ(truth.size(), 7u);
  expect_all_found(truth, run(truth), 1);
}

TEST(Pipeline, HandlesMeasurementSpaceEdges) {
  const auto truth = load_targets_csv(std::string(RP_SCENARIO_DIR) + "/edge_cases.csv");
  expect_all_found(truth, run(truth), 1);
}

TEST(Pipeline, NoiseOnlyFrameProducesFewFalseAlarms) {
  std::size_t total = 0;
  for (std::uint32_t seed = 1; seed <= 5; ++seed) total += run({}, seed).detections.size();
  EXPECT_LE(total, 2u);
}

TEST(Pipeline, TracksGroundTruthOverFrameSequence) {
  // Targets move between frames; every frame must still match ground truth.
  RadarConfig cfg;
  Simulator sim(cfg);
  RadarPipeline pipe(cfg);
  const std::vector<Target> initial{{30.0, 15.0, 0.1, 10.0}, {180.0, -30.0, -0.2, 20.0}};
  DataCube adc;
  FrameResult res;
  for (int f = 0; f < 10; ++f) {
    std::vector<Target> truth = initial;
    for (auto& t : truth) t.range_m += t.velocity_mps * 0.1 * f;
    sim.simulate_frame(truth, adc);
    pipe.process(adc, res);
    const auto ev = evaluate(truth, res.detections, {1.0, 1.0});
    EXPECT_EQ(ev.false_negatives, 0u) << "frame " << f;
  }
}

TEST(Pipeline, IsDeterministicForFixedSeed) {
  const std::vector<Target> truth{{50.0, -10.0, 0.2, 10.0}, {120.0, 5.0, -0.3, 15.0}};
  const auto a = run(truth, 9);
  const auto b = run(truth, 9);
  ASSERT_EQ(a.detections.size(), b.detections.size());
  for (std::size_t i = 0; i < a.detections.size(); ++i) {
    EXPECT_EQ(a.detections[i].range_m, b.detections[i].range_m);
    EXPECT_EQ(a.detections[i].velocity_mps, b.detections[i].velocity_mps);
    EXPECT_EQ(a.detections[i].azimuth_rad, b.detections[i].azimuth_rad);
  }
}

TEST(Pipeline, ReportsStageTimings) {
  const auto res = run({{50.0, 0.0, 0.0, 10.0}});
  EXPECT_GT(res.timings.total_ms, 0.0);
  EXPECT_GE(res.timings.total_ms, res.timings.range_doppler_ms + res.timings.cfar_ms);
}

}  // namespace radar
