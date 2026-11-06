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
#include <random>
#include <sstream>
#include <string>

#include "radar/calibration.hpp"
#include "radar/clustering.hpp"
#include "radar/ego_estimator.hpp"
#include "radar/ego_velocity.hpp"
#include "radar/evaluation.hpp"
#include "radar/fusion.hpp"
#include "radar/impairment_experiments.hpp"
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

// TDM-MIMO 3 Tx x 4 Rx: random targets across the full unfolded velocity span
// (+/-70 m/s, three times the TDM Doppler span) and +/-12 deg elevation.
TEST(Regression, Mimo) {
  const RadarConfig radar = RadarConfig::tdm_mimo_3tx4rx();
  SimulatorConfig sc;
  sc.seed = 21;
  Simulator sim(radar, sc);
  RadarPipeline pipe(radar);
  PipelineConfig no_comp;
  no_comp.tdm_motion_compensation = false;
  RadarPipeline pipe_uncompensated(radar, no_comp);

  std::mt19937 rng(5);
  std::uniform_real_distribution<double> vel(-65.0, 65.0), az(-45.0, 45.0), el(-12.0, 12.0),
      rcs(5.0, 15.0);
  std::size_t tp = 0, fn = 0, fp = 0, unfolded = 0, moving = 0;
  double er = 0, ev = 0, ea = 0, ee = 0, ea_uncomp = 0;
  std::size_t tp_uncomp = 0;
  DataCube adc;
  FrameResult res;
  const int frames = 40;
  for (int f = 0; f < frames; ++f) {
    // Four targets, 25 m apart in range so they never share a cell.
    std::vector<Target> truth;
    for (int i = 0; i < 4; ++i) {
      Target t;
      t.range_m = 15.0 + 25.0 * i + 10.0 * std::uniform_real_distribution<double>(0, 1)(rng);
      t.velocity_mps = vel(rng);
      t.azimuth_rad = az(rng) * kDeg;
      t.elevation_rad = el(rng) * kDeg;
      t.rcs_dbsm = rcs(rng);
      truth.push_back(t);
    }
    sim.simulate_frame(truth, adc);
    pipe.process(adc, res);
    MatchGates gates;
    gates.velocity_mps = 1.0;  // a wrong unfold is 46 m/s off: counts as missed
    const auto e = evaluate(truth, res.detections, gates);
    tp += e.true_positives;
    fn += e.false_negatives;
    fp += e.false_positives;
    const double n = static_cast<double>(e.true_positives);
    er += e.mean_abs_range_error_m * n;
    ev += e.mean_abs_velocity_error_mps * n;
    ea += e.mean_abs_azimuth_error_deg * n;
    ee += e.mean_abs_elevation_error_deg * n;
    for (std::size_t i = 0; i < truth.size(); ++i) {
      if (std::abs(truth[i].velocity_mps) < radar.max_velocity()) continue;
      ++moving;  // only these need unfolding
      if (e.truth_to_detection[i] >= 0) ++unfolded;
    }

    pipe_uncompensated.process(adc, res);
    MatchGates loose;
    loose.velocity_mps = 2.0 * radar.max_unfolded_velocity();  // velocity is folded here
    const auto u = evaluate(truth, res.detections, loose);
    tp_uncomp += u.true_positives;
    ea_uncomp += u.mean_abs_azimuth_error_deg * static_cast<double>(u.true_positives);
  }
  const double n = static_cast<double>(tp);
  check("mimo.detection_rate", n / static_cast<double>(tp + fn));
  check("mimo.false_alarms_per_frame", static_cast<double>(fp) / frames);
  check("mimo.unfold_success_rate", static_cast<double>(unfolded) / static_cast<double>(moving));
  check("mimo.range_mae_m", er / n);
  check("mimo.velocity_mae_mps", ev / n);
  check("mimo.azimuth_mae_deg", ea / n);
  check("mimo.elevation_mae_deg", ee / n);
  check("mimo.uncompensated_azimuth_mae_deg", ea_uncomp / static_cast<double>(tp_uncomp));
}

// Impairments and mitigations (fixed seeds; see impairments_demo for the
// full-length runs and figures).
TEST(Regression, Interference) {
  const auto r = experiments::run_interference(12, 11);
  check("interference.floor_rise_db", r.floor_rise_db);
  check("interference.detection_rate_unmitigated", r.interfered.detection_rate);
  check("interference.detection_rate_mitigated", r.mitigated.detection_rate);
  check("interference.false_alarms_mitigated", r.mitigated.false_alarms_per_frame);
  check("interference.residual_floor_db", r.residual_floor_db);
  check("interference.weak_target_rate_mitigated", r.weak_target_rate_mitigated);
}

TEST(Regression, ClutterCfar) {
  const auto r = experiments::run_clutter(10, 11);
  check("clutter.ca_masked_rate", r.ca.masked_rate);
  check("clutter.ca_clutter_edge_rate", r.ca.clutter_edge_rate);
  check("clutter.os_masked_rate", r.os.masked_rate);
  check("clutter.os_clutter_edge_rate", r.os.clutter_edge_rate);
  check("clutter.os_false_alarms_per_frame", r.os.false_alarms_per_frame);
}

TEST(Regression, Multipath) {
  const auto r = experiments::run_multipath(20, 11);
  check("multipath.mirror_ghost_removal",
        1.0 - r.mirror_ghosts_left_per_frame / r.mirror_ghosts_per_frame);
  check("multipath.real_flagged_per_frame", r.real_flagged_per_frame);
  check("multipath.rail_offset_error_m", r.rail_offset_error_m);
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
  double vel_sq = 0.0, fused_sq = 0.0, imu_sq = 0.0;
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
      const Vec2 p = sf.ego.pose.translation();
      fused_sq += (fused.pose().translation() - p).squaredNorm();
      imu_sq += (imu_only.pose().translation() - p).squaredNorm();
    }
    if (k < steps) scene.step(0.01);
  }
  // RMSE over the run rather than the final error: a single end-point sample
  // swings with tiny floating-point differences between compilers.
  const Pose2D truth = scene.ego().pose;
  const double fused_rmse = std::sqrt(fused_sq / static_cast<double>(frames));
  const double imu_rmse = std::sqrt(imu_sq / static_cast<double>(frames));
  check("odometry.ego_velocity_valid_fraction",
        static_cast<double>(valid) / static_cast<double>(frames));
  check("odometry.ego_velocity_rmse_mps", std::sqrt(vel_sq / static_cast<double>(vel_n)));
  check("odometry.fused_position_rmse_m", fused_rmse);
  check("odometry.fused_final_yaw_error_deg",
        std::abs(wrap_angle(fused.pose().yaw - truth.yaw)) / kDeg);
  check("odometry.imu_only_to_fused_rmse_ratio", imu_rmse / std::max(1e-3, fused_rmse));
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
