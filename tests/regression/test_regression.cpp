// Scenario-level regression suite.
//
// Every experiment runs with fixed seeds and reports metrics that are checked
// against tests/regression/thresholds.txt ("<metric> <min|max> <value>").
// A Markdown report of measured values vs. bounds is written to
// regression_report.md in the working directory, for CI artifacts.
//
// Set RP_SKIP_TIMING=1 to skip latency budgets (emulators, sanitizers).

#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "radar/calibration.hpp"
#include "radar/clustering.hpp"
#include "radar/ego_estimator.hpp"
#include "radar/ego_velocity.hpp"
#include "radar/evaluation.hpp"
#include "radar/fusion.hpp"
#include "radar/imu.hpp"
#include "radar/io.hpp"
#include "radar/pipeline.hpp"
#include "radar/scenarios.hpp"
#include "radar/simulator.hpp"
#include "radar/tracker.hpp"
#include "radar/tracking_metrics.hpp"

namespace radar {
namespace {

constexpr double kDeg = kPi / 180.0;

struct Bound {
  bool is_min = true;
  double value = 0.0;
};

struct ReportRow {
  std::string metric;
  double measured;
  Bound bound;
  bool pass;
};

std::map<std::string, Bound>& thresholds() {
  static std::map<std::string, Bound> t = [] {
    std::map<std::string, Bound> m;
    std::ifstream f(RP_THRESHOLDS_FILE);
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string name, kind;
      double v;
      if (ss >> name >> kind >> v) m[name] = {kind == "min", v};
    }
    return m;
  }();
  return t;
}

std::vector<ReportRow>& report() {
  static std::vector<ReportRow> r;
  return r;
}

void check(const std::string& metric, double measured) {
  std::printf("[ metric   ] %-42s %.6g\n", metric.c_str(), measured);
  const auto it = thresholds().find(metric);
  ASSERT_NE(it, thresholds().end()) << "no threshold for " << metric;
  const Bound b = it->second;
  const bool pass = b.is_min ? measured >= b.value : measured <= b.value;
  report().push_back({metric, measured, b, pass});
  EXPECT_TRUE(pass) << metric << " = " << measured << (b.is_min ? " < min " : " > max ") << b.value;
}

bool skip_timing() {
  const char* v = std::getenv("RP_SKIP_TIMING");
  return v && std::string(v) != "0";
}

class ReportWriter : public ::testing::Environment {
 public:
  void TearDown() override {
    std::ofstream f("regression_report.md");
    f << "| Metric | Measured | Bound | |\n|---|---|---|---|\n";
    for (const auto& r : report()) {
      f << "| `" << r.metric << "` | " << r.measured << " | " << (r.bound.is_min ? "≥ " : "≤ ")
        << r.bound.value << " | " << (r.pass ? "pass" : "**FAIL**") << " |\n";
    }
  }
};
const auto* const kReport = ::testing::AddGlobalTestEnvironment(new ReportWriter);

double calibrated_rcs_offset(const RadarConfig& radar, RadarPipeline& pipeline) {
  SimulatorConfig sc;
  sc.seed = 999;
  Simulator sim(radar, sc);
  FrameResult res;
  double sum = 0.0;
  for (int i = 0; i < 5; ++i) {
    pipeline.process(sim.simulate_frame({{50.0, 0.0, 0.0, 10.0}}), res);
    double best = -1e9;
    for (const auto& d : res.detections) best = std::max(best, static_cast<double>(d.snr_db));
    sum += best;
  }
  return calibrate_rcs_offset(sum / 5.0, 50.0, 10.0);
}

}  // namespace

TEST(Regression, PointCloud) {
  RadarConfig radar;
  Simulator sim(radar);
  RadarPipeline pipe(radar);
  const auto scenario = load_targets_csv(std::string(RP_SCENARIO_DIR) + "/multi_target.csv");
  std::size_t tp = 0, fn = 0, fp = 0;
  double er = 0, ev = 0, ea = 0, latency = 0;
  DataCube adc;
  FrameResult res;
  const int frames = 20;
  for (int f = 0; f < frames; ++f) {
    std::vector<Target> truth;
    for (Target t : scenario) {
      t.range_m += t.velocity_mps * 0.05 * f;
      truth.push_back(t);
    }
    sim.simulate_frame(truth, adc);
    pipe.process(adc, res);
    latency += res.timings.total_ms;
    const auto e = evaluate(truth, res.detections);
    tp += e.true_positives;
    fn += e.false_negatives;
    fp += e.false_positives;
    er += e.mean_abs_range_error_m * static_cast<double>(e.true_positives);
    ev += e.mean_abs_velocity_error_mps * static_cast<double>(e.true_positives);
    ea += e.mean_abs_azimuth_error_deg * static_cast<double>(e.true_positives);
  }
  const double n = static_cast<double>(tp);
  check("point_cloud.detection_rate", n / static_cast<double>(tp + fn));
  check("point_cloud.false_alarms_per_frame", static_cast<double>(fp) / frames);
  check("point_cloud.range_mae_m", er / n);
  check("point_cloud.velocity_mae_mps", ev / n);
  check("point_cloud.azimuth_mae_deg", ea / n);
  if (!skip_timing()) check("point_cloud.latency_mean_ms", latency / frames);
}

TEST(Regression, Tracking) {
  RadarConfig radar;
  SimulatorConfig sc;
  sc.seed = 7;
  Simulator sim(radar, sc);
  RadarPipeline pipe(radar);
  pipe.set_rcs_offset_db(calibrated_rcs_offset(radar, pipe));
  Scene scene = scenarios::intersection();
  Dbscan db;
  MultiTargetTracker trk;
  TrackingEvaluator eval;
  DataCube adc;
  FrameResult res;
  std::vector<Cluster> cl;
  for (int f = 0; f < 120; ++f) {  // 6 s
    const SceneFrame sf = scene.snapshot();
    sim.simulate_frame(sf.targets, adc);
    pipe.process(adc, res);
    db.cluster(res.detections, cl);
    trk.update(sf.time_s, cl, {sf.sensor_pose, sf.sensor_velocity});
    eval.add_frame(sf.objects, trk.tracks());
    scene.step(0.05);
  }
  const auto& m = eval.metrics();
  const auto veh = static_cast<Eigen::Index>(ObjectClass::kVehicle);
  const double veh_total = static_cast<double>(m.confusion.row(veh).sum());
  check("tracking.mota", m.mota());
  check("tracking.recall", m.recall());
  check("tracking.id_switches", static_cast<double>(m.id_switches));
  check("tracking.velocity_error_mps", m.mean_velocity_error_mps());
  check("tracking.vehicle_class_accuracy",
        static_cast<double>(m.confusion(veh, veh)) / std::max(1.0, veh_total));
}

TEST(Regression, Odometry) {
  const double duration = 15.0;
  Scene scene = scenarios::ego_drive({}, duration);
  const Pose2D mount = scene.config().sensor_mount;
  RadarConfig radar;
  SimulatorConfig sc;
  sc.seed = 3;
  Simulator sim(radar, sc);
  RadarPipeline pipe(radar);
  ImuConfig ic;
  ic.seed = 13;
  ImuSimulator imu(ic);
  EgoEstimatorConfig ec;
  ec.sensor_mount = mount;
  EgoEstimator fused(ec), imu_only(ec);
  EgoVelocityEstimator ego_vel;
  fused.initialize(scene.ego().pose, scene.ego().velocity(), 0.0);
  imu_only.initialize(scene.ego().pose, scene.ego().velocity(), 0.0);

  DataCube adc;
  FrameResult res;
  double vel_sq = 0.0;
  std::size_t vel_n = 0, valid = 0, frames = 0;
  const auto steps = static_cast<int>(duration / 0.01);
  for (int k = 0; k <= steps; ++k) {
    const ImuSample s = imu.sample(scene.ego(), 0.01 * k);
    fused.propagate(s);
    imu_only.propagate(s);
    if (k % 5 == 0) {
      const SceneFrame sf = scene.snapshot();
      sim.simulate_frame(sf.targets, adc);
      pipe.process(adc, res);
      ++frames;
      const auto e = ego_vel.estimate(res.detections);
      if (e.valid) {
        ++valid;
        const Vec2 truth = sf.sensor_pose.rotation_matrix().transpose() * sf.sensor_velocity;
        vel_sq += (e.velocity - truth).squaredNorm();
        ++vel_n;
        fused.update_radar_velocity(e.velocity, e.covariance);
      }
    }
    if (k < steps) scene.step(0.01);
  }
  const Pose2D truth = scene.ego().pose;
  const double fused_err = (fused.pose().translation() - truth.translation()).norm();
  const double imu_err = (imu_only.pose().translation() - truth.translation()).norm();
  check("odometry.ego_velocity_valid_fraction",
        static_cast<double>(valid) / static_cast<double>(frames));
  check("odometry.ego_velocity_rmse_mps", std::sqrt(vel_sq / static_cast<double>(vel_n)));
  check("odometry.fused_final_error_m", fused_err);
  check("odometry.fused_final_yaw_error_deg",
        std::abs(wrap_angle(fused.pose().yaw - truth.yaw)) / kDeg);
  check("odometry.imu_only_to_fused_error_ratio", imu_err / std::max(1e-3, fused_err));
  check("odometry.gyro_bias_error_radps", std::abs(fused.gyro_bias() - imu.gyro_bias()));
}

TEST(Regression, ExtrinsicCalibration) {
  CameraModel truth;
  truth.T_cam_from_radar =
      make_camera_extrinsic(Vec3(-1.8, 0.1, 0.8), 0.4 * kDeg, 2.0 * kDeg, -0.7 * kDeg);
  const Pose3 rough = make_camera_extrinsic(Vec3(-1.65, 0.0, 0.9), -1.0 * kDeg, 0.0, 1.8 * kDeg);
  RadarConfig radar;
  SimulatorConfig sc;
  sc.seed = 82;
  Simulator sim(radar, sc);
  RadarPipeline pipe(radar);
  std::mt19937 rng(5);
  std::normal_distribution<double> px(0.0, 1.0);
  std::vector<Correspondence> data;
  FrameResult res;
  int i = 0;
  for (double range : {8.0, 12.0, 18.0, 25.0}) {
    for (double az_deg : {-25.0, -10.0, 5.0, 20.0}) {
      const double az = az_deg * kDeg, z = 0.3 + 0.35 * (i++ % 4) - 0.6;
      pipe.process(sim.simulate_frame({{range, 0.0, az, 15.0}}), res);
      ASSERT_FALSE(res.detections.empty());
      const auto& d = *std::max_element(
          res.detections.begin(), res.detections.end(),
          [](const Detection& a, const Detection& b) { return a.snr_db < b.snr_db; });
      Correspondence c;
      c.p_radar = Vec3(d.x_m, d.y_m, z);
      ASSERT_TRUE(truth.project(Vec3(range * std::cos(az), range * std::sin(az), z), c.uv));
      c.uv += Eigen::Vector2d(px(rng), px(rng));
      data.push_back(c);
    }
  }
  const auto cal = calibrate_extrinsic(truth.K, data, rough);
  const PoseError e = pose_difference(truth.T_cam_from_radar, cal.T_cam_from_radar);
  check("calibration.rotation_error_deg", e.rotation_rad / kDeg);
  check("calibration.translation_error_m", e.translation_m);
  check("calibration.reprojection_rms_px", cal.final_rms_px);
}

}  // namespace radar
