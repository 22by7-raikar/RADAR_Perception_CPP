// End-to-end demo: simulate a sequence of frames, run the perception pipeline
// on each, score the point cloud against ground truth and report latency.
//
//   radar_demo [--scenario FILE] [--frames N] [--out DIR] [--seed S]

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "radar/evaluation.hpp"
#include "radar/io.hpp"
#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"

namespace {

struct Args {
  std::string scenario;
  std::string out_dir = "out";
  std::size_t frames = 50;
  double frame_period_s = 0.05;  // 20 Hz
  std::uint32_t seed = 42;
};

void usage() {
  std::puts("usage: radar_demo [--scenario FILE] [--frames N] [--out DIR] [--seed S]");
}

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
    else if (k == "--frames")
      a.frames = std::stoul(next());
    else if (k == "--out")
      a.out_dir = next();
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::stoul(next()));
    else if (k == "-h" || k == "--help") {
      usage();
      std::exit(0);
    } else
      throw std::invalid_argument("unknown argument " + k);
  }
  return a;
}

std::vector<radar::Target> default_scenario() {
  constexpr double deg = radar::kPi / 180.0;
  return {
      {12.0, 0.0, -25.0 * deg, 15.0},   // static structure, close
      {35.0, -4.5, 10.0 * deg, 0.0},    // slow, small, approaching
      {60.0, 12.0, -8.0 * deg, 20.0},   // receding vehicle
      {60.0, -12.0, 15.0 * deg, 20.0},  // same range, opposite velocity
      {95.0, -25.0, 30.0 * deg, 10.0}, {150.0, 8.0, -40.0 * deg, 25.0},
      {210.0, -40.0, 5.0 * deg, 20.0},  // far and fast
  };
}

std::vector<radar::Target> propagate(const std::vector<radar::Target>& initial, double t,
                                     double max_range) {
  std::vector<radar::Target> out;
  for (radar::Target tgt : initial) {
    tgt.range_m += tgt.velocity_mps * t;
    if (tgt.range_m > 1.0 && tgt.range_m < max_range) out.push_back(tgt);
  }
  return out;
}

double percentile(std::vector<double> v, double p) {
  if (v.empty()) return 0.0;
  std::sort(v.begin(), v.end());
  const auto idx = static_cast<std::size_t>(p * static_cast<double>(v.size() - 1) + 0.5);
  return v[idx];
}

double mean(const std::vector<double>& v) {
  return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Args args = parse(argc, argv);

    radar::RadarConfig radar;
    radar::SimulatorConfig sim_cfg;
    sim_cfg.seed = args.seed;
    radar::PipelineConfig pipe_cfg;

    const auto scenario =
        args.scenario.empty() ? default_scenario() : radar::load_targets_csv(args.scenario);

    std::printf("FMCW radar: %.0f GHz, %zu samples x %zu chirps x %zu rx\n", radar.carrier_hz / 1e9,
                radar.samples_per_chirp, radar.chirps_per_frame, radar.num_rx);
    std::printf("  range res %.2f m, max range %.0f m | velocity res %.2f m/s, max +/-%.1f m/s\n",
                radar.range_resolution(), radar.max_range(), radar.velocity_resolution(),
                radar.max_velocity());
    std::printf("  CFAR: Pfa %.0e, train %zux%zu, guard %zux%zu\n\n", pipe_cfg.cfar.pfa,
                pipe_cfg.cfar.train_range, pipe_cfg.cfar.train_doppler, pipe_cfg.cfar.guard_range,
                pipe_cfg.cfar.guard_doppler);

    radar::Simulator sim(radar, sim_cfg);
    radar::RadarPipeline pipeline(radar, pipe_cfg);
    std::filesystem::create_directories(args.out_dir);
    const std::string det_path = args.out_dir + "/detections.csv";
    const std::string truth_path = args.out_dir + "/truth.csv";

    std::ofstream timing_csv(args.out_dir + "/timings.csv");
    timing_csv
        << "frame,range_doppler_ms,integration_ms,cfar_ms,peak_grouping_ms,angle_ms,total_ms\n";

    radar::DataCube frame;
    radar::FrameResult result;
    std::size_t tp = 0, fn = 0, fp = 0;
    double err_r = 0.0, err_v = 0.0, err_a = 0.0;
    std::vector<double> t_total, t_rd, t_int, t_cfar, t_peak, t_ang;

    for (std::size_t f = 0; f < args.frames; ++f) {
      const double t = static_cast<double>(f) * args.frame_period_s;
      const auto truth = propagate(scenario, t, radar.max_range());
      sim.simulate_frame(truth, frame);
      pipeline.process(frame, result);

      const auto ev = radar::evaluate(truth, result.detections);
      tp += ev.true_positives;
      fn += ev.false_negatives;
      fp += ev.false_positives;
      const double n = static_cast<double>(ev.true_positives);
      err_r += ev.mean_abs_range_error_m * n;
      err_v += ev.mean_abs_velocity_error_mps * n;
      err_a += ev.mean_abs_azimuth_error_deg * n;

      t_total.push_back(result.timings.total_ms);
      t_rd.push_back(result.timings.range_doppler_ms);
      t_int.push_back(result.timings.integration_ms);
      t_cfar.push_back(result.timings.cfar_ms);
      t_peak.push_back(result.timings.peak_grouping_ms);
      t_ang.push_back(result.timings.angle_ms);

      const auto& tm = result.timings;
      timing_csv << f << ',' << tm.range_doppler_ms << ',' << tm.integration_ms << ',' << tm.cfar_ms
                 << ',' << tm.peak_grouping_ms << ',' << tm.angle_ms << ',' << tm.total_ms << '\n';

      radar::write_detections_csv(det_path, result.detections, f, f > 0);
      radar::write_targets_csv(truth_path, truth, f, f > 0);

      if (f == 0) {
        radar::write_power_map_csv(args.out_dir + "/power_map_db.csv", pipeline.power_map());
        std::printf("Frame 0: %zu CFAR cells -> %zu detections\n", result.cfar_cells,
                    result.detections.size());
        std::printf("  %8s %9s %8s %8s %8s %7s\n", "range_m", "vel_mps", "az_deg", "x_m", "y_m",
                    "snr_db");
        for (const auto& d : result.detections) {
          std::printf("  %8.2f %9.2f %8.2f %8.2f %8.2f %7.1f\n", d.range_m, d.velocity_mps,
                      d.azimuth_rad * 180.0 / radar::kPi, d.x_m, d.y_m, d.snr_db);
        }
        std::puts("");
      }
    }

    const double tpd = std::max<double>(1.0, static_cast<double>(tp));
    std::printf("Over %zu frames:\n", args.frames);
    std::printf(
        "  detection rate     %.1f %% (%zu / %zu)\n",
        100.0 * static_cast<double>(tp) / std::max<double>(1.0, static_cast<double>(tp + fn)), tp,
        tp + fn);
    std::printf(
        "  false alarms/frame %.2f\n",
        static_cast<double>(fp) / static_cast<double>(std::max<std::size_t>(1, args.frames)));
    std::printf("  mean |error|       range %.3f m, velocity %.3f m/s, azimuth %.2f deg\n",
                err_r / tpd, err_v / tpd, err_a / tpd);
    std::printf("\nLatency per frame (ms): mean %.3f  p50 %.3f  p99 %.3f\n", mean(t_total),
                percentile(t_total, 0.5), percentile(t_total, 0.99));
    std::printf(
        "  range-Doppler FFT %.3f | integration %.3f | CFAR %.3f | peaks %.3f | angle %.3f\n",
        mean(t_rd), mean(t_int), mean(t_cfar), mean(t_peak), mean(t_ang));
    std::printf("\nOutputs written to %s/\n", args.out_dir.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    usage();
    return 1;
  }
}
