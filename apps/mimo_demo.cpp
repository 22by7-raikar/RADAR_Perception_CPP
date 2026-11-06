// TDM-MIMO demo (3 Tx x 4 Rx, 12 virtual channels).
//
// 1. Velocity sweep: one target at 25 deg azimuth, range rate -68 .. +68 m/s.
//    Azimuth error and reported velocity with no TDM motion compensation,
//    compensation without unfolding, and compensation + unfolding.
// 2. Overhead structure: at 100 km/h a sign gantry over the road, a stopped
//    car under it and a metal plate on the road all show up as static
//    reflectors in the lane. A radar without elevation has to treat all of
//    them as obstacles (or ignore all static returns); with elevation only
//    the car is in the drivable height band.
//
//   mimo_demo [--out DIR] [--seed N]

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"

namespace {

using namespace radar;

constexpr double kDeg = kPi / 180.0;

struct Args {
  std::string out_dir = "out_mimo";
  std::uint32_t seed = 3;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::invalid_argument("missing value for " + k);
      return argv[++i];
    };
    if (k == "--out")
      a.out_dir = next();
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::stoul(next()));
    else
      throw std::invalid_argument("unknown argument " + k);
  }
  return a;
}

/// Strongest detection within 3 m of the expected range, or nullptr.
const Detection* strongest_near(const std::vector<Detection>& dets, double range) {
  const Detection* best = nullptr;
  for (const auto& d : dets) {
    if (std::abs(d.range_m - range) > 3.0) continue;
    if (!best || d.snr_db > best->snr_db) best = &d;
  }
  return best;
}

void velocity_sweep(const Args& args) {
  const RadarConfig radar = RadarConfig::tdm_mimo_3tx4rx();
  PipelineConfig none, comp_only, full;
  none.tdm_motion_compensation = false;
  comp_only.velocity_unfolding = false;
  RadarPipeline p_none(radar, none), p_comp(radar, comp_only), p_full(radar, full);
  SimulatorConfig sc;
  sc.seed = args.seed;
  Simulator sim(radar, sc);

  std::ofstream f(args.out_dir + "/velocity_sweep.csv");
  f << "velocity,az_err_none,az_err_comp,az_err_full,v_comp,v_full\n";
  const double az = 25.0;
  std::size_t full_ok = 0, comp_ok = 0, total = 0;
  DataCube adc;
  FrameResult r_none, r_comp, r_full;
  for (double v = -68.0; v <= 68.0 + 1e-9; v += 1.0) {
    Target t;
    t.range_m = 60.0;
    t.velocity_mps = v;
    t.azimuth_rad = az * kDeg;
    sim.simulate_frame({t}, adc);
    p_none.process(adc, r_none);
    p_comp.process(adc, r_comp);
    p_full.process(adc, r_full);
    const Detection* a = strongest_near(r_none.detections, t.range_m);
    const Detection* b = strongest_near(r_comp.detections, t.range_m);
    const Detection* c = strongest_near(r_full.detections, t.range_m);
    auto err = [&](const Detection* d) {
      return d ? std::abs(d->azimuth_rad / kDeg - az) : std::nan("");
    };
    auto vel = [](const Detection* d) { return d ? d->velocity_mps : std::nanf(""); };
    f << v << ',' << err(a) << ',' << err(b) << ',' << err(c) << ',' << vel(b) << ',' << vel(c)
      << '\n';
    ++total;
    if (c && std::abs(c->velocity_mps - v) < 0.5 && err(c) < 1.0) ++full_ok;
    if (b && std::abs(b->velocity_mps - v) < 0.5 && err(b) < 1.0) ++comp_ok;
  }
  std::printf("Velocity sweep (-68 .. 68 m/s, TDM Doppler span +/-%.1f m/s):\n",
              radar.max_velocity());
  std::printf("  correct velocity and azimuth: compensation only %zu/%zu, + unfolding %zu/%zu\n",
              comp_ok, total, full_ok, total);
}

struct SceneObject {
  std::string name;
  double x, y, height;  ///< relative to the radar footprint on the ground [m]
  double rcs_dbsm;
  double speed;  ///< along +x, world frame [m/s]
};

/// Decision for a static in-lane return: the vehicle can pass under it, pass
/// over it, or must treat it as an obstacle. Only call pass-under/over when
/// the height is outside the vehicle's band by two standard deviations;
/// without elevation the uncertainty is infinite and every in-lane return is
/// an obstacle.
constexpr double kLaneHalfWidth = 1.8;
constexpr double kDriveOverHeight = 0.25;  ///< below: plates, joints, cat's eyes
constexpr double kClearanceHeight = 2.5;   ///< above: signs, gantries, bridges
constexpr double kRadarHeight = 0.7;
constexpr double kEgoSpeed = 27.8;  // 100 km/h

bool in_path_obstacle(const Detection& d) {
  if (std::abs(d.y_m) >= kLaneHalfWidth) return false;
  const double h = static_cast<double>(d.z_m) + kRadarHeight;
  const double h_std = static_cast<double>(d.range_m) * static_cast<double>(d.elevation_std_rad);
  if (h - 2.0 * h_std > kClearanceHeight) return false;  // pass under
  if (h + 2.0 * h_std < kDriveOverHeight) return false;  // drive over
  return true;
}

void overhead_scene(const Args& args) {
  std::vector<SceneObject> objects;
  for (double y : {-6.0, -3.0, 0.0, 3.0, 6.0}) objects.push_back({"gantry", 70.0, y, 5.8, 8.0, 0});
  objects.push_back({"car", 105.0, 0.2, 0.6, 12.0, 0.0});
  objects.push_back({"car", 105.6, -0.6, 0.9, 6.0, 0.0});
  objects.push_back({"plate", 42.0, -0.4, 0.0, 2.0, 0.0});
  objects.push_back({"lead", 55.0, -3.6, 0.6, 10.0, 22.0});

  struct Variant {
    std::string name;
    RadarConfig radar;
  };
  const std::vector<Variant> variants{{"3d", RadarConfig{}},
                                      {"4d", RadarConfig::tdm_mimo_3tx4rx()}};

  std::ofstream f(args.out_dir + "/overhead_detections.csv");
  f << "radar,frame,object,x,y,height,height_std,velocity,obstacle\n";
  std::printf("\nOverhead structure at %.0f km/h (static world closes at %.1f m/s):\n",
              kEgoSpeed * 3.6, kEgoSpeed);
  std::printf("  %-6s %-14s %-14s %-14s %-14s\n", "radar", "gantry", "stopped car", "road plate",
              "adjacent lead");
  for (const Variant& var : variants) {
    SimulatorConfig sc;
    sc.seed = args.seed;
    Simulator sim(var.radar, sc);
    RadarPipeline pipe(var.radar);
    DataCube adc;
    FrameResult res;
    constexpr int kFrames = 10;
    std::size_t flagged_gantry = 0, flagged_car = 0, flagged_plate = 0, flagged_lead = 0;
    for (int fr = 0; fr < kFrames; ++fr) {
      const double t = 0.05 * fr;
      std::vector<Target> targets;
      std::vector<std::string> labels;
      std::vector<double> ranges;
      for (const auto& o : objects) {
        const double x = o.x + (o.speed - kEgoSpeed) * t;
        const double dz = o.height - kRadarHeight;
        const double ground = std::hypot(x, o.y);
        Target tg;
        tg.range_m = std::hypot(ground, dz);
        tg.azimuth_rad = std::atan2(o.y, x);
        tg.elevation_rad = std::atan2(dz, ground);
        tg.velocity_mps = (o.speed - kEgoSpeed) * x / tg.range_m;
        tg.rcs_dbsm = o.rcs_dbsm;
        targets.push_back(tg);
        labels.push_back(o.name);
        ranges.push_back(tg.range_m);
      }
      sim.simulate_frame(targets, adc);
      pipe.process(adc, res);
      bool g = false, c = false, p = false, l = false;
      for (const auto& d : res.detections) {
        // Label each detection with the nearest object in range.
        std::size_t best = 0;
        for (std::size_t i = 1; i < ranges.size(); ++i) {
          if (std::abs(ranges[i] - d.range_m) < std::abs(ranges[best] - d.range_m)) best = i;
        }
        const std::string& name = labels[best];
        const bool obstacle = in_path_obstacle(d);
        f << var.name << ',' << fr << ',' << name << ',' << d.x_m << ',' << d.y_m << ','
          << d.z_m + kRadarHeight << ',' << d.range_m * d.elevation_std_rad << ',' << d.velocity_mps
          << ',' << obstacle << '\n';
        if (!obstacle) continue;
        g |= name == "gantry";
        c |= name == "car";
        p |= name == "plate";
        l |= name == "lead";
      }
      flagged_gantry += g;
      flagged_car += c;
      flagged_plate += p;
      flagged_lead += l;
    }
    auto cell = [](std::size_t n) {
      return std::to_string(n) + "/" + std::to_string(kFrames) + " frames";
    };
    std::printf("  %-6s %-14s %-14s %-14s %-14s\n", var.name == "3d" ? "3D" : "4D",
                cell(flagged_gantry).c_str(), cell(flagged_car).c_str(),
                cell(flagged_plate).c_str(), cell(flagged_lead).c_str());
  }
  std::printf("  (frames in which the object raised an in-lane obstacle detection)\n");
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Args args = parse(argc, argv);
    std::filesystem::create_directories(args.out_dir);
    velocity_sweep(args);
    overhead_scene(args);
    std::printf("\nCSV written to %s/\n", args.out_dir.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
