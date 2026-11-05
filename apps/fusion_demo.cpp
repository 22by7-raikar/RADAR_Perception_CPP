// Radar-camera fusion demo.
//
// 1. Extrinsic calibration: corner reflectors at known heights are measured by
//    the radar (full signal chain) and the camera; LM on SE(3) recovers the
//    radar-to-camera transform from a rough mounting guess.
// 2. Drive: radar tracks (20 Hz) are fused with camera detections (25 Hz,
//    unsynchronised). Four configurations are compared on the same data:
//    calibrated + time sync, calibrated without time sync, rough extrinsic,
//    and radar only.
//
//   fusion_demo [--duration S] [--out DIR] [--seed N]

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "radar/calibration.hpp"
#include "radar/clustering.hpp"
#include "radar/fusion.hpp"
#include "radar/pipeline.hpp"
#include "radar/scenarios.hpp"
#include "radar/simulator.hpp"
#include "radar/tracker.hpp"

namespace {

using namespace radar;

struct Args {
  double duration_s = 20.0;
  std::string out_dir = "out_fusion";
  std::uint32_t seed = 5;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + k);
      return argv[++i];
    };
    if (k == "--duration")
      a.duration_s = std::stod(next());
    else if (k == "--out")
      a.out_dir = next();
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::stoul(next()));
    else
      throw std::invalid_argument("unknown argument " + k);
  }
  return a;
}

constexpr double kDeg = kPi / 180.0;
constexpr double kRadarHeight = 0.6;

/// Calibration session: one reflector per placement, measured through the radar chain.
std::vector<Correspondence> collect_calibration(const RadarConfig& radar, RadarPipeline& pipeline,
                                                const CameraModel& truth_camera,
                                                std::uint32_t seed) {
  SimulatorConfig sc;
  sc.seed = seed + 77;
  Simulator sim(radar, sc);
  std::mt19937 rng(seed);
  std::normal_distribution<double> px(0.0, 1.0);
  std::vector<Correspondence> out;
  FrameResult res;
  for (double range : {8.0, 12.0, 18.0, 25.0}) {
    for (double az_deg : {-25.0, -10.0, 5.0, 20.0}) {
      const double az = az_deg * kDeg;
      const double height = 0.3 + 0.35 * static_cast<double>(out.size() % 4);  // tripod setting
      pipeline.process(sim.simulate_frame({{range, 0.0, az, 15.0}}), res);
      if (res.detections.empty()) continue;
      const auto best = *std::max_element(
          res.detections.begin(), res.detections.end(),
          [](const Detection& a, const Detection& b) { return a.snr_db < b.snr_db; });
      Correspondence c;
      c.p_radar = Vec3(best.x_m, best.y_m, height - kRadarHeight);
      const Vec3 truth(range * std::cos(az), range * std::sin(az), height - kRadarHeight);
      if (!truth_camera.project(truth, c.uv)) continue;
      c.uv += Eigen::Vector2d(px(rng), px(rng));
      out.push_back(c);
    }
  }
  return out;
}

struct ConfigStats {
  std::size_t objects = 0;  ///< fused objects whose ground truth is visible to the camera
  std::size_t correct = 0;  ///< associated to the camera detection of the same object
  std::size_t wrong = 0;    ///< associated to another detection
  double iou_sum = 0.0;
  std::size_t class_correct[2] = {0, 0};  ///< pedestrian, vehicle
  std::size_t class_total[2] = {0, 0};

  void print(const char* name) const {
    const double n = static_cast<double>(std::max<std::size_t>(1, objects));
    auto pct = [](std::size_t a, std::size_t b) { return b ? 100.0 * a / b : 0.0; };
    std::printf("  %-34s %6.1f %%  %5.1f %%  %5.2f   %5.1f %%   %5.1f %%\n", name,
                100.0 * correct / n, 100.0 * wrong / n,
                correct ? iou_sum / static_cast<double>(correct) : 0.0,
                pct(class_correct[0], class_total[0]), pct(class_correct[1], class_total[1]));
  }
};

/// Tracker output and sensor state at one radar frame, kept for latency experiments.
struct RadarSnapshot {
  double time = 0.0;
  std::vector<Track> tracks;
  SensorState sensor;
};

}  // namespace

int main(int argc, char** argv) {
  try {
    const Args args = parse(argc, argv);
    RadarConfig radar;
    SimulatorConfig sc;
    sc.seed = args.seed;
    Simulator sim(radar, sc);
    RadarPipeline pipeline(radar);

    // --- Extrinsic calibration -------------------------------------------
    CameraModel truth_cam;
    truth_cam.T_cam_from_radar =
        make_camera_extrinsic(Vec3(-1.8, 0.1, 0.8), 0.4 * kDeg, 2.0 * kDeg, -0.7 * kDeg);
    CameraModel rough_cam = truth_cam;  // tape-measure / CAD guess
    rough_cam.T_cam_from_radar =
        make_camera_extrinsic(Vec3(-1.65, 0.0, 0.9), -1.0 * kDeg, 0.0 * kDeg, 1.8 * kDeg);

    const auto corr = collect_calibration(radar, pipeline, truth_cam, args.seed);
    const CalibrationResult cal =
        calibrate_extrinsic(truth_cam.K, corr, rough_cam.T_cam_from_radar);
    CameraModel calib_cam = truth_cam;
    calib_cam.T_cam_from_radar = cal.T_cam_from_radar;
    const PoseError e_rough =
        pose_difference(truth_cam.T_cam_from_radar, rough_cam.T_cam_from_radar);
    const PoseError e_cal = pose_difference(truth_cam.T_cam_from_radar, cal.T_cam_from_radar);
    std::printf("Extrinsic calibration from %zu reflector placements (%d LM iterations)\n",
                corr.size(), cal.iterations);
    std::printf("  reprojection RMS  %.1f px -> %.2f px\n", cal.initial_rms_px, cal.final_rms_px);
    std::printf("  rotation error    %.2f deg -> %.3f deg\n", e_rough.rotation_rad / kDeg,
                e_cal.rotation_rad / kDeg);
    std::printf("  position error    %.3f m -> %.3f m\n\n", e_rough.translation_m,
                e_cal.translation_m);

    // --- Drive -------------------------------------------------------------
    SceneConfig scfg;
    Scene scene = scenarios::ego_drive(scfg, args.duration_s);
    CameraSimulator camera(truth_cam, kRadarHeight);
    FusionConfig synced;
    FusionConfig unsynced;
    unsynced.compensate_time = false;
    // Fusion is stateful (per-track class evidence), so each variant owns one.
    std::vector<RadarCameraFusion> fusions = {
        {calib_cam, synced},   {rough_cam, synced}, {calib_cam, synced},  {calib_cam, synced},
        {calib_cam, unsynced}, {calib_cam, synced}, {calib_cam, unsynced}};

    // A variant fuses camera frame t with the radar snapshot from t - latency.
    // With time sync, tracks are predicted to t and placed with the sensor pose
    // at t; without it, the stale tracks are used as they are, with the stale pose.
    struct Variant {
      const char* name;
      RadarCameraFusion* fusion;
      double latency_s;
      bool time_sync;
      bool radar_only;
      ConfigStats stats;
    };
    std::vector<Variant> variants = {
        {"calibrated + time sync", &fusions[0], 0.0, true, false, {}},
        {"rough extrinsic + time sync", &fusions[1], 0.0, true, false, {}},
        {"radar only (no camera)", &fusions[2], 0.0, true, true, {}},
        {"100 ms latency, time sync", &fusions[3], 0.1, true, false, {}},
        {"100 ms latency, no time sync", &fusions[4], 0.1, false, false, {}},
        {"250 ms latency, time sync", &fusions[5], 0.25, true, false, {}},
        {"250 ms latency, no time sync", &fusions[6], 0.25, false, false, {}},
    };

    Dbscan dbscan;
    MultiTargetTracker tracker;
    DataCube adc;
    FrameResult res;
    std::vector<Cluster> clusters;
    std::vector<RadarSnapshot> history;

    std::filesystem::create_directories(args.out_dir);
    std::ofstream f_cam(args.out_dir + "/camera_frame.csv");
    f_cam << "kind,id,class,u0,v0,u1,v1,label,matched\n";
    bool frame_written = false;

    const double dt = 0.005;
    const auto steps = static_cast<int>(args.duration_s / dt);
    for (int k = 0; k <= steps; ++k) {
      const double t = k * dt;
      if (k % 10 == 0) {  // radar, 20 Hz
        const SceneFrame sf = scene.snapshot();
        sim.simulate_frame(sf.targets, adc);
        pipeline.process(adc, res);
        dbscan.cluster(res.detections, clusters);
        const SensorState s{sf.sensor_pose, sf.sensor_velocity};
        tracker.update(t, clusters, s);
        history.push_back({t, tracker.tracks(), s});
        if (history.size() > 20) history.erase(history.begin());
      }
      if (k % 8 == 4 && t > 0.5) {  // camera, 25 Hz, offset from the radar clock
        const SceneFrame cf = scene.snapshot();
        const auto dets = camera.detect(cf);
        const SensorState s_cam{cf.sensor_pose, cf.sensor_velocity};

        // Ground truth for a fused object: nearest camera-visible object,
        // compared at the camera time.
        auto truth_of = [&](const Vec2& p_sensor) -> const GroundTruthObject* {
          const GroundTruthObject* best = nullptr;
          double bd = 1e9;
          for (const auto& g : cf.objects) {
            if (g.cls != ObjectClass::kVehicle && g.cls != ObjectClass::kPedestrian) continue;
            const double d = (g.position_sensor - p_sensor).norm();
            const double gate = 1.5 + 0.5 * std::max(g.length, g.width);
            if (d < gate && d < bd) {
              bd = d;
              best = &g;
            }
          }
          if (!best) return nullptr;
          for (const auto& d : dets) {
            if (d.truth_id == best->id) return best;  // only score objects the camera saw
          }
          return nullptr;
        };

        for (Variant& v : variants) {
          const RadarSnapshot* snap = nullptr;
          for (const auto& h : history) {
            if (h.time <= t - v.latency_s + 1e-9) snap = &h;
          }
          if (!snap) continue;
          const SensorState& pose = v.time_sync ? s_cam : snap->sensor;
          const auto fused = v.fusion->fuse(snap->tracks, snap->time, pose, t, dets);
          for (const auto& f : fused) {
            // Score against where the object really is at camera time.
            const Vec2 p_now = s_cam.pose.to_local(pose.pose.to_parent(f.position_sensor));
            const GroundTruthObject* g = truth_of(v.time_sync ? f.position_sensor : p_now);
            if (!g) continue;
            ++v.stats.objects;
            ObjectClass predicted = f.cls;
            if (v.radar_only) {
              for (const auto& tr : snap->tracks) {
                if (tr.id == f.track_id) predicted = tr.cls;
              }
            } else if (f.camera_index >= 0) {
              const auto& d = dets[static_cast<std::size_t>(f.camera_index)];
              if (d.truth_id == g->id) {
                ++v.stats.correct;
                v.stats.iou_sum += f.iou;
              } else {
                ++v.stats.wrong;
              }
            }
            const int ci = g->cls == ObjectClass::kPedestrian ? 0 : 1;
            ++v.stats.class_total[ci];
            if (predicted == g->cls) ++v.stats.class_correct[ci];
          }

          if (&v == &variants.front() && !frame_written && t > 6.5 && fused.size() >= 3) {
            for (const auto& d : dets) {
              f_cam << "camera," << d.truth_id << ',' << to_string(d.cls) << ',' << d.box.u0 << ','
                    << d.box.v0 << ',' << d.box.u1 << ',' << d.box.v1 << ",\n";
            }
            for (const auto& f : fused) {
              char label[64];
              std::snprintf(label, sizeof(label), "%.0f m, %.1f m/s", f.position_sensor.norm(),
                            f.velocity_world.norm());
              f_cam << "track," << f.track_id << ',' << to_string(f.cls) << ',' << f.projected.u0
                    << ',' << f.projected.v0 << ',' << f.projected.u1 << ',' << f.projected.v1
                    << ',' << label << ',' << (f.camera_index >= 0) << '\n';
            }
            for (const auto& d : res.detections) {
              Eigen::Vector2d uv;
              if (calib_cam.project(Vec3(d.x_m, d.y_m, 0.0), uv)) {
                f_cam << "point,-1,," << uv.x() << ',' << uv.y() << ',' << uv.x() << ',' << uv.y()
                      << ',' << d.range_m << '\n';
              }
            }
            frame_written = true;
          }
        }
      }
      scene.step(dt);
    }

    std::printf("Fusion over %.0f s (radar 20 Hz, camera 25 Hz, unsynchronised clocks)\n",
                args.duration_s);
    std::printf("  %-34s %8s %8s %6s %9s %9s\n", "configuration", "correct", "wrong", "IoU",
                "ped cls", "veh cls");
    for (const Variant& v : variants) v.stats.print(v.name);
    std::printf("Outputs written to %s/\n", args.out_dir.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
