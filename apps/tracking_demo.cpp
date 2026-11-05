// Object-level demo: scene -> IQ frames -> point cloud -> clusters -> tracks,
// scored against ground truth with CLEAR-MOT metrics and a class confusion matrix.
//
//   tracking_demo [--scenario intersection|ego_drive] [--duration S] [--out DIR] [--seed N]

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "radar/clustering.hpp"
#include "radar/pipeline.hpp"
#include "radar/scenarios.hpp"
#include "radar/simulator.hpp"
#include "radar/tracker.hpp"
#include "radar/tracking_metrics.hpp"

namespace {

struct Args {
  std::string scenario = "intersection";
  double duration_s = 10.0;
  double frame_period_s = 0.05;
  std::string out_dir = "out_tracking";
  std::uint32_t seed = 7;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + k);
      return argv[++i];
    };
    if (k == "--scenario")
      a.scenario = next();
    else if (k == "--duration")
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

}  // namespace

namespace radar {

/// Measures the RCS offset with a simulated 10 dBsm reflector at 50 m, the way
/// a corner reflector is used to calibrate a real sensor.
double calibrate(const RadarConfig& radar, RadarPipeline& pipeline, std::uint32_t seed) {
  SimulatorConfig sc;
  sc.seed = seed + 1000;
  Simulator sim(radar, sc);
  FrameResult res;
  double snr_sum = 0.0;
  constexpr int kShots = 10;
  for (int i = 0; i < kShots; ++i) {
    pipeline.process(sim.simulate_frame({{50.0, 0.0, 0.0, 10.0}}), res);
    double best = -1e9;
    for (const auto& d : res.detections) best = std::max(best, static_cast<double>(d.snr_db));
    snr_sum += best;
  }
  return calibrate_rcs_offset(snr_sum / kShots, 50.0, 10.0);
}

}  // namespace radar

int main(int argc, char** argv) {
  using namespace radar;
  try {
    const Args args = parse(argc, argv);
    Scene scene = args.scenario == "ego_drive" ? scenarios::ego_drive() : scenarios::intersection();

    RadarConfig radar;
    SimulatorConfig sc;
    sc.seed = args.seed;
    Simulator sim(radar, sc);
    RadarPipeline pipeline(radar);
    const double offset = calibrate(radar, pipeline, args.seed);
    pipeline.set_rcs_offset_db(offset);
    std::printf("RCS calibration offset: %.2f dB\n", offset);

    Dbscan dbscan;
    MultiTargetTracker tracker;
    TrackingEvaluator evaluator;

    std::filesystem::create_directories(args.out_dir);
    std::ofstream f_trk(args.out_dir + "/tracks.csv");
    std::ofstream f_gt(args.out_dir + "/truth_objects.csv");
    std::ofstream f_det(args.out_dir + "/detections_world.csv");
    std::ofstream f_ego(args.out_dir + "/ego.csv");
    f_trk << "frame,time,id,confirmed,x,y,vx,vy,class,confidence\n";
    f_gt << "frame,time,id,class,x,y,vx,vy\n";
    f_det << "frame,x,y,velocity,cluster,rcs_dbsm\n";
    f_ego << "frame,time,x,y,yaw,speed\n";

    DataCube adc;
    FrameResult res;
    std::vector<Cluster> clusters;
    double t_cluster = 0.0, t_track = 0.0, t_signal = 0.0;
    const auto frames = static_cast<std::size_t>(args.duration_s / args.frame_period_s);

    for (std::size_t f = 0; f < frames; ++f) {
      const SceneFrame sf = scene.snapshot();
      sim.simulate_frame(sf.targets, adc);
      pipeline.process(adc, res);
      t_signal += res.timings.total_ms;

      const auto c0 = std::chrono::steady_clock::now();
      dbscan.cluster(res.detections, clusters);
      const auto c1 = std::chrono::steady_clock::now();
      tracker.update(sf.time_s, clusters, {sf.sensor_pose, sf.sensor_velocity});
      const auto c2 = std::chrono::steady_clock::now();
      t_cluster += std::chrono::duration<double, std::milli>(c1 - c0).count();
      t_track += std::chrono::duration<double, std::milli>(c2 - c1).count();

      evaluator.add_frame(sf.objects, tracker.tracks());

      for (const auto& t : tracker.tracks()) {
        f_trk << f << ',' << sf.time_s << ',' << t.id << ',' << t.confirmed() << ',' << t.x(0)
              << ',' << t.x(1) << ',' << t.x(2) << ',' << t.x(3) << ',' << to_string(t.cls) << ','
              << t.class_confidence << '\n';
      }
      for (const auto& g : sf.objects) {
        f_gt << f << ',' << sf.time_s << ',' << g.id << ',' << to_string(g.cls) << ','
             << g.position_world.x() << ',' << g.position_world.y() << ',' << g.velocity_world.x()
             << ',' << g.velocity_world.y() << '\n';
      }
      const auto& labels = dbscan.labels();
      for (std::size_t i = 0; i < res.detections.size(); ++i) {
        const auto& d = res.detections[i];
        const Vec2 w = sf.sensor_pose.to_parent({d.x_m, d.y_m});
        f_det << f << ',' << w.x() << ',' << w.y() << ',' << d.velocity_mps << ',' << labels[i]
              << ',' << d.rcs_dbsm << '\n';
      }
      f_ego << f << ',' << sf.time_s << ',' << sf.ego.pose.x << ',' << sf.ego.pose.y << ','
            << sf.ego.pose.yaw << ',' << sf.ego.speed << '\n';

      scene.step(args.frame_period_s);
    }

    const TrackingMetrics& m = evaluator.metrics();
    const double n = static_cast<double>(frames);
    std::printf("\nScenario '%s', %zu frames (%.1f s)\n", args.scenario.c_str(), frames,
                args.duration_s);
    std::printf("  MOTA %.3f | recall %.3f | MOTP %.2f m | velocity error %.2f m/s\n", m.mota(),
                m.recall(), m.motp_m(), m.mean_velocity_error_mps());
    std::printf("  matches %zu / %zu, misses %zu, false tracks %zu, ID switches %zu\n", m.matches,
                m.ground_truth, m.misses, m.false_positives, m.id_switches);
    std::printf("  class accuracy %.1f %%\n", 100.0 * m.class_accuracy());
    std::printf(
        "  confusion (rows = truth, cols = predicted: unknown pedestrian vehicle static)\n");
    for (int r = 1; r <= 2; ++r) {
      std::printf("    %-10s", to_string(static_cast<ObjectClass>(r)));
      for (int c = 0; c < kNumObjectClasses; ++c) std::printf(" %6zu", m.confusion(r, c));
      std::printf("\n");
    }
    std::printf("  mean latency per frame: signal %.2f ms, clustering %.3f ms, tracking %.3f ms\n",
                t_signal / n, t_cluster / n, t_track / n);
    std::printf("Outputs written to %s/\n", args.out_dir.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
