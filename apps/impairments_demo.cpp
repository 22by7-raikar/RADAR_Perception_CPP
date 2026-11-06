// Real-world impairments and their mitigations.
//
// 1. Interference: another FMCW radar's chirps sweep through the IF band.
//    Clean vs. interfered vs. time-domain blanking.
// 2. Clutter and masking: vegetation and a strong target next to weak ones.
//    Cell-averaging vs. ordered-statistic CFAR.
// 3. Multipath: guardrail reflections create ghost targets. Guardrail fitted
//    from static detections; mirrored ghosts behind it are flagged.
//
//   impairments_demo [--frames N] [--out DIR] [--seed N]

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "radar/impairment_experiments.hpp"

namespace {

using namespace radar;
using namespace radar::experiments;

struct Args {
  int frames = 40;
  std::string out_dir = "out_impairments";
  std::uint32_t seed = 11;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + k);
      return argv[++i];
    };
    if (k == "--frames")
      a.frames = std::stoi(next());
    else if (k == "--out")
      a.out_dir = next();
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::stoul(next()));
    else
      throw std::invalid_argument("unknown argument " + k);
  }
  return a;
}

void interference(const Args& a) {
  const InterferenceResult r = run_interference(a.frames, a.seed);
  std::printf("1. Interference (INR 30 dB per sample, %d frames)\n", a.frames);
  std::printf("   %-26s %10s %14s %16s\n", "", "detected", "false alarms", "weak targets");
  std::printf("   %-26s %9.1f%% %14.2f %15.1f%%\n", "clean", 100 * r.clean.detection_rate,
              r.clean.false_alarms_per_frame, 100 * r.weak_target_rate_clean);
  std::printf("   %-26s %9.1f%% %14.2f %15.1f%%\n", "interfered", 100 * r.interfered.detection_rate,
              r.interfered.false_alarms_per_frame, 100 * r.weak_target_rate_interfered);
  std::printf("   %-26s %9.1f%% %14.2f %15.1f%%\n", "interfered + blanking",
              100 * r.mitigated.detection_rate, r.mitigated.false_alarms_per_frame,
              100 * r.weak_target_rate_mitigated);
  std::printf("   noise floor +%.1f dB; after blanking %+.2f dB; %.2f%% of samples blanked\n",
              r.floor_rise_db, r.residual_floor_db, 100 * r.blanked_fraction);
  std::printf("   frame time %.2f ms -> %.2f ms with mitigation\n\n", r.frame_ms_plain,
              r.frame_ms_mitigated);

  std::ofstream f(a.out_dir + "/interference_chirp.csv");
  f << "sample,envelope,residual,blanked\n";
  for (std::size_t n = 0; n < r.chirp_envelope.size(); ++n)
    f << n << ',' << r.chirp_envelope[n] << ',' << r.chirp_residual[n] << ','
      << int(r.chirp_blanked[n]) << '\n';
  std::ofstream g(a.out_dir + "/interference_profile.csv");
  g << "range_bin,clean_db,interfered_db,mitigated_db\n";
  for (std::size_t n = 0; n < r.range_profile_clean_db.size(); ++n)
    g << n << ',' << r.range_profile_clean_db[n] << ',' << r.range_profile_interfered_db[n] << ','
      << r.range_profile_mitigated_db[n] << '\n';
}

void clutter(const Args& a) {
  const ClutterResult r = run_clutter(a.frames, a.seed);
  std::printf("2. Clutter and masking (%d frames)\n", a.frames);
  std::printf("   %-8s %12s %14s %11s %16s %13s %9s\n", "CFAR", "clutter edge", "behind strong",
              "reference", "clutter det/frm", "false alarms", "CFAR ms");
  for (const auto& [name, c] : {std::pair<const char*, CfarComparison>{"CA", r.ca}, {"OS", r.os}}) {
    std::printf("   %-8s %11.0f%% %13.0f%% %10.0f%% %16.1f %13.2f %9.2f\n", name,
                100 * c.clutter_edge_rate, 100 * c.masked_rate, 100 * c.reference_rate,
                c.clutter_detections_per_frame, c.false_alarms_per_frame, c.cfar_ms);
  }
  std::printf("\n");
  std::ofstream f(a.out_dir + "/clutter_power.csv");
  f << "rows," << r.rows << ",cols," << r.cols << '\n';
  for (std::size_t i = 0; i < r.power_db.size(); ++i)
    f << r.power_db[i] << ((i + 1) % r.cols == 0 ? '\n' : ',');
  std::ofstream t(a.out_dir + "/clutter_truth.csv");
  t << "name,range,velocity\n";
  const char* names[] = {"pedestrian", "car", "motorbike", "reference"};
  for (std::size_t i = 0; i < r.truth.size(); ++i)
    t << names[i] << ',' << r.truth[i].range_m << ',' << r.truth[i].velocity_mps << '\n';
  std::ofstream g(a.out_dir + "/clutter_detections.csv");
  g << "cfar,range,velocity,azimuth,snr\n";
  for (const auto& [name, dets] :
       {std::pair<const char*, const std::vector<Detection>*>{"CA", &r.ca_detections},
        {"OS", &r.os_detections}}) {
    for (const Detection& d : *dets)
      g << name << ',' << d.range_m << ',' << d.velocity_mps << ',' << d.azimuth_rad << ','
        << d.snr_db << '\n';
  }
}

void multipath(const Args& a) {
  const MultipathResult r = run_multipath(a.frames, a.seed);
  std::printf("3. Multipath off a guardrail (%d frames)\n", a.frames);
  std::printf("   guardrail found in %.0f%% of frames, offset error %.2f m\n",
              100 * r.rail_found_fraction, r.rail_offset_error_m);
  std::printf("   mirror ghosts per frame: %.2f -> %.2f after the filter\n",
              r.mirror_ghosts_per_frame, r.mirror_ghosts_left_per_frame);
  std::printf("   real moving detections wrongly flagged: %.2f of %.2f per frame\n",
              r.real_flagged_per_frame, r.real_moving_per_frame);
  std::printf("   direct-bearing single-bounce ghosts (not addressed): %.2f per frame\n\n",
              r.direct_ghosts_per_frame);
  std::ofstream f(a.out_dir + "/multipath_detections.csv");
  f << "x,y,velocity,source,flagged\n";
  const char* names[] = {"real", "rail", "ghost_mirror", "ghost_direct"};
  for (std::size_t i = 0; i < r.detections.size(); ++i) {
    f << r.detections[i].x_m << ',' << r.detections[i].y_m << ',' << r.detections[i].velocity_mps
      << ',' << names[static_cast<int>(r.sources[i])] << ',' << int(r.flagged[i]) << '\n';
  }
  std::ofstream g(a.out_dir + "/multipath_rail.csv");
  g << "valid,slope,offset\n"
    << r.rail.valid << ',' << r.rail.slope << ',' << r.rail.offset << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Args args = parse(argc, argv);
    std::filesystem::create_directories(args.out_dir);
    interference(args);
    clutter(args);
    multipath(args);
    std::printf("CSV written to %s/\n", args.out_dir.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
