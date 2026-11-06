#include "radar/impairment_experiments.hpp"

#include <algorithm>
#include <cmath>
#include <random>

#include "radar/evaluation.hpp"
#include "radar/impairments.hpp"
#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"

namespace radar::experiments {

namespace {

constexpr double kDeg = kPi / 180.0;
constexpr double kFrameDt = 0.05;

Target target(double range, double velocity, double az_deg, double rcs) {
  Target t;
  t.range_m = range;
  t.velocity_mps = velocity;
  t.azimuth_rad = az_deg * kDeg;
  t.rcs_dbsm = rcs;
  return t;
}

std::vector<float> mean_over_doppler_db(const Grid2D<float>& power) {
  std::vector<float> out(power.rows());
  for (std::size_t r = 0; r < power.rows(); ++r) {
    double s = 0.0;
    for (std::size_t d = 0; d < power.cols(); ++d) s += power(r, d);
    out[r] = static_cast<float>(10.0 * std::log10(s / static_cast<double>(power.cols())));
  }
  return out;
}

}  // namespace

// --- 1. Interference ---------------------------------------------------------

InterferenceResult run_interference(int frames, std::uint32_t seed) {
  const RadarConfig radar;
  SimulatorConfig sc;
  sc.seed = seed;
  Simulator sim(radar, sc);
  PipelineConfig on;
  on.interference.enabled = true;
  RadarPipeline plain(radar), mitigated(radar, on);

  // The multi-target scene plus three weak (0 dBsm) targets beyond 100 m.
  const std::vector<Target> base{target(12, 0.0, -25, 15),  target(35, -4.5, 10, 0),
                                 target(60, 12.0, -8, 20),  target(60, -12.0, 15, 20),
                                 target(95, -25.0, 30, 10), target(150, 8.0, -40, 25),
                                 target(210, -40.0, 5, 20), target(120, -3.0, -15, 0),
                                 target(150, 5.0, 20, 0),   target(135, -10.0, 40, 0)};
  auto is_weak = [](const Target& t) { return t.rcs_dbsm <= 0.0 && t.range_m > 100.0; };

  std::mt19937 rng(seed + 100);
  std::uniform_real_distribution<double> offset(0.0, 37e-6), az(-40.0, 40.0);

  InterferenceResult res;
  std::size_t tp[3] = {0, 0, 0}, fp[3] = {0, 0, 0}, weak_hit[3] = {0, 0, 0};
  double ms_plain = 0.0, ms_mitigated = 0.0;
  std::size_t total = 0, weak_total = 0, blanked = 0;
  double rise = 0.0, residual = 0.0;
  DataCube clean, dirty;
  FrameResult rc, ri, rm;
  for (int f = 0; f < frames; ++f) {
    std::vector<Target> truth;
    for (Target t : base) {
      t.range_m += t.velocity_mps * kFrameDt * f;
      truth.push_back(t);
    }
    sim.simulate_frame(truth, clean);
    dirty = clean;
    Interferer in;
    in.time_offset_s = offset(rng);
    in.azimuth_rad = az(rng) * kDeg;
    add_interference(radar, in, dirty);

    plain.process(clean, rc);
    const double floor_clean = rc.noise_floor;
    const auto ec = evaluate(truth, rc.detections);
    plain.process(dirty, ri);
    ms_plain += ri.timings.total_ms;
    const auto ei = evaluate(truth, ri.detections);
    rise += 10.0 * std::log10(ri.noise_floor / floor_clean);
    if (f == frames - 1) res.range_profile_interfered_db = mean_over_doppler_db(plain.power_map());
    mitigated.process(dirty, rm);
    ms_mitigated += rm.timings.total_ms;
    const auto em = evaluate(truth, rm.detections);
    residual += 10.0 * std::log10(rm.noise_floor / floor_clean);
    blanked += rm.interference_samples_blanked;
    if (f == frames - 1) {
      res.range_profile_mitigated_db = mean_over_doppler_db(mitigated.power_map());
      plain.process(clean, rc);
      res.range_profile_clean_db = mean_over_doppler_db(plain.power_map());
    }

    const EvaluationResult* e[3] = {&ec, &ei, &em};
    for (int k = 0; k < 3; ++k) {
      tp[k] += e[k]->true_positives;
      fp[k] += e[k]->false_positives;
    }
    total += truth.size();
    for (std::size_t i = 0; i < truth.size(); ++i) {
      if (!is_weak(truth[i])) continue;
      ++weak_total;
      weak_hit[0] += ec.truth_to_detection[i] >= 0;
      weak_hit[1] += ei.truth_to_detection[i] >= 0;
      weak_hit[2] += em.truth_to_detection[i] >= 0;
    }

    if (f == frames - 1) {
      // A chirp the burst hits, for the time-domain figure.
      for (std::size_t m = 0; m < radar.chirps_per_frame; ++m) {
        std::vector<float> env(radar.samples_per_chirp, 0.0f);
        bool hit = false;
        for (std::size_t r = 0; r < radar.num_rx; ++r) {
          for (std::size_t n = 0; n < radar.samples_per_chirp; ++n) {
            env[n] += squared_magnitude(dirty(r, m, n));
            hit |= dirty(r, m, n) != clean(r, m, n);
          }
        }
        if (!hit) continue;
        res.chirp_envelope = env;
        std::vector<float> scratch;
        DataCube one(radar.num_rx, 1, radar.samples_per_chirp);
        for (std::size_t r = 0; r < radar.num_rx; ++r)
          for (std::size_t n = 0; n < radar.samples_per_chirp; ++n) one(r, 0, n) = dirty(r, m, n);
        InterferenceMitigationConfig zeroing = on.interference;
        zeroing.reconstruct_iterations = 0;  // keep the blanked slots visible
        mitigate_interference(radar, zeroing, one, scratch);
        res.chirp_residual =
            notched_residual(&dirty(0, m, 0), radar.samples_per_chirp, on.interference.notch_db);
        res.chirp_blanked.resize(radar.samples_per_chirp);
        for (std::size_t n = 0; n < radar.samples_per_chirp; ++n)
          res.chirp_blanked[n] = one(0, 0, n) == cfloat{};
        break;
      }
    }
  }
  const auto n = static_cast<double>(total);
  DetectionStats* s[3] = {&res.clean, &res.interfered, &res.mitigated};
  for (int k = 0; k < 3; ++k) {
    s[k]->detection_rate = static_cast<double>(tp[k]) / n;
    s[k]->false_alarms_per_frame = static_cast<double>(fp[k]) / frames;
  }
  res.floor_rise_db = rise / frames;
  res.residual_floor_db = residual / frames;
  res.blanked_fraction =
      static_cast<double>(blanked) /
      (static_cast<double>(frames) * static_cast<double>(radar.chirps_per_frame) *
       static_cast<double>(radar.samples_per_chirp));
  const auto weak_n = static_cast<double>(weak_total);
  res.weak_target_rate_clean = static_cast<double>(weak_hit[0]) / weak_n;
  res.weak_target_rate_interfered = static_cast<double>(weak_hit[1]) / weak_n;
  res.weak_target_rate_mitigated = static_cast<double>(weak_hit[2]) / weak_n;
  res.frame_ms_plain = ms_plain / frames;
  res.frame_ms_mitigated = ms_mitigated / frames;
  return res;
}

// --- 2. Clutter --------------------------------------------------------------

ClutterResult run_clutter(int frames, std::uint32_t seed) {
  const RadarConfig radar;  // static sensor
  SimulatorConfig sc;
  sc.seed = seed;
  Simulator sim(radar, sc);
  PipelineConfig ca_cfg, os_cfg;
  os_cfg.cfar_type = CfarType::kOrderedStatistic;
  RadarPipeline ca(radar, ca_cfg), os(radar, os_cfg);

  // Roadside vegetation 70-110 m, 25-50 deg to the right, moving in the wind.
  ClutterPatch patch;
  patch.range_min_m = 70.0;
  patch.range_max_m = 110.0;
  patch.azimuth_min_rad = -50.0 * kDeg;
  patch.azimuth_max_rad = -25.0 * kDeg;
  patch.rcs_dbsm = -6.0;
  const std::vector<Target> clutter = clutter_targets(patch, seed + 7);

  ClutterResult res;
  struct Counts {
    std::size_t edge = 0, masked = 0, reference = 0, clutter = 0, fa = 0;
    double ms = 0.0;
  } counts[2];
  DataCube adc;
  FrameResult r;
  for (int f = 0; f < frames; ++f) {
    const double t = kFrameDt * f;
    // Pedestrian (-8 dBsm) walking out of the vegetation towards the road,
    // just in front of it in range; a car with a motorbike 4 m behind it at a
    // similar speed; the same motorbike-sized target alone as the reference.
    const Target edge = target(67.0 - 1.3 * t, -1.3, -35.0, -8.0);
    const Target strong = target(100.0 - 9.0 * t, -9.0, 10.0, 25.0);
    const Target masked = target(104.0 - 9.3 * t, -9.3, 14.0, -1.0);
    const Target reference = target(104.0 - 9.3 * t, 9.3, -10.0, -1.0);
    std::vector<Target> scene = clutter;
    scene.insert(scene.end(), {edge, strong, masked, reference});
    sim.simulate_frame(scene, adc);

    RadarPipeline* pipes[2] = {&ca, &os};
    std::vector<Detection> frame_dets[2];
    bool os_found_masked = false;
    for (int k = 0; k < 2; ++k) {
      pipes[k]->process(adc, r);
      Counts& c = counts[k];
      c.ms += r.timings.cfar_ms;
      auto near = [](const Detection& d, const Target& tg) {
        return std::abs(d.range_m - tg.range_m) < 1.5 &&
               std::abs(d.velocity_mps - tg.velocity_mps) < 1.0;
      };
      bool e = false, m = false, ref = false;
      for (const Detection& d : r.detections) {
        if (near(d, edge)) {
          e = true;
        } else if (near(d, masked)) {
          m = true;
        } else if (near(d, reference)) {
          ref = true;
        } else if (near(d, strong)) {
        } else if (d.range_m > patch.range_min_m - 1.5 && d.range_m < patch.range_max_m + 1.5 &&
                   std::abs(d.velocity_mps) < 1.5) {
          ++c.clutter;
        } else {
          ++c.fa;
        }
      }
      c.edge += e;
      c.masked += m;
      c.reference += ref;
      frame_dets[k] = r.detections;
      if (k == 1) os_found_masked = m;
    }
    // Figure snapshot: the last frame that shows the typical OS outcome.
    if (os_found_masked || res.power_db.empty()) {
      res.ca_detections = frame_dets[0];
      res.os_detections = frame_dets[1];
      res.truth = {edge, strong, masked, reference};
      const auto& p = os.power_map();
      res.rows = p.rows();
      res.cols = p.cols();
      res.power_db.resize(p.size());
      for (std::size_t i = 0; i < p.size(); ++i)
        res.power_db[i] = 10.0f * std::log10(std::max(p.data()[i], 1e-12f));
    }
  }
  CfarComparison* out[2] = {&res.ca, &res.os};
  for (int k = 0; k < 2; ++k) {
    const auto n = static_cast<double>(frames);
    out[k]->clutter_edge_rate = static_cast<double>(counts[k].edge) / n;
    out[k]->masked_rate = static_cast<double>(counts[k].masked) / n;
    out[k]->reference_rate = static_cast<double>(counts[k].reference) / n;
    out[k]->clutter_detections_per_frame = static_cast<double>(counts[k].clutter) / n;
    out[k]->false_alarms_per_frame = static_cast<double>(counts[k].fa) / n;
    out[k]->cfar_ms = counts[k].ms / n;
  }
  return res;
}

// --- 3. Multipath ------------------------------------------------------------

MultipathResult run_multipath(int frames, std::uint32_t seed) {
  const RadarConfig radar;
  SimulatorConfig sc;
  sc.seed = seed;
  Simulator sim(radar, sc);
  RadarPipeline pipe(radar);
  const double ego_speed = 15.0;
  const Vec2 ego_v(ego_speed, 0.0);
  PlanarReflector rail;
  rail.y_m = -3.5;
  rail.loss_db = 3.0;
  GhostFilterConfig gcfg;

  MultipathResult res;
  double mirror = 0, mirror_left = 0, direct = 0, real_flagged = 0, real_moving = 0, rail_err = 0;
  int rail_found = 0;
  DataCube adc;
  FrameResult r;
  for (int f = 0; f < frames; ++f) {
    const double t = kFrameDt * f;
    std::vector<Target> targets;
    std::vector<Source> source;
    // Guardrail posts every 2 m, static in the world.
    for (double x = 6.0 - std::fmod(ego_speed * t, 2.0); x < 150.0; x += 2.0) {
      targets.push_back(direct_target({{x, rail.y_m}, {-ego_speed, 0.0}, -6.0}));
      source.push_back(Source::kRail);
    }
    // Moving objects (sensor frame, velocity relative to the sensor).
    const std::vector<MovingPoint> objects{
        {{35.0 - 4.0 * t, 0.0}, {-4.0, 0.0}, 12.0},     // lead car, slower
        {{70.0 + 5.0 * t, -1.0}, {5.0, 0.0}, 10.0},     // car pulling away
        {{110.0 - 35.0 * t, 3.5}, {-35.0, 0.0}, 18.0},  // oncoming truck
        {{60.0 + 2.0 * t, -9.0}, {2.0, 0.0}, 3.0},      // motorbike behind the rail
    };
    for (const auto& o : objects) {
      targets.push_back(direct_target(o));
      source.push_back(Source::kReal);
      const auto mp = multipath_targets(o, rail);
      for (std::size_t i = 0; i < mp.size(); ++i) {
        targets.push_back(mp[i]);
        source.push_back(i == 1 ? Source::kGhostDirect : Source::kGhostMirror);
      }
    }
    sim.simulate_frame(targets, adc);
    pipe.process(adc, r);

    // Attribute each detection to the nearest simulated path.
    std::vector<Source> label(r.detections.size(), Source::kRail);
    for (std::size_t i = 0; i < r.detections.size(); ++i) {
      const Detection& d = r.detections[i];
      double best = 1e18;
      for (std::size_t j = 0; j < targets.size(); ++j) {
        const double dr = (d.range_m - targets[j].range_m) / 1.0;
        const double dv = (d.velocity_mps - targets[j].velocity_mps) / 0.5;
        const double da = (d.azimuth_rad - targets[j].azimuth_rad) / (3.0 * kDeg);
        const double cost = dr * dr + dv * dv + da * da;
        if (cost < best) {
          best = cost;
          label[i] = source[j];
        }
      }
    }

    std::vector<ReflectorLine> lines{fit_reflector(r.detections, ego_v, -1, gcfg),
                                     fit_reflector(r.detections, ego_v, +1, gcfg)};
    std::vector<char> ghost;
    flag_multipath_ghosts(r.detections, lines, ego_v, gcfg, ghost);
    if (lines[0].valid) {
      ++rail_found;
      rail_err += std::abs(lines[0].offset - rail.y_m);
    }
    for (std::size_t i = 0; i < label.size(); ++i) {
      switch (label[i]) {
        case Source::kGhostMirror:
          ++mirror;
          mirror_left += !ghost[i];
          break;
        case Source::kGhostDirect:
          ++direct;
          break;
        case Source::kReal:
          ++real_moving;
          real_flagged += ghost[i];
          break;
        case Source::kRail:
          break;
      }
    }
    if (f == frames - 1) {
      res.detections = r.detections;
      res.sources = label;
      res.flagged = ghost;
      res.rail = lines[0];
    }
  }
  const auto n = static_cast<double>(frames);
  res.mirror_ghosts_per_frame = mirror / n;
  res.mirror_ghosts_left_per_frame = mirror_left / n;
  res.direct_ghosts_per_frame = direct / n;
  res.real_flagged_per_frame = real_flagged / n;
  res.real_moving_per_frame = real_moving / n;
  res.rail_found_fraction = rail_found / n;
  res.rail_offset_error_m = rail_found > 0 ? rail_err / rail_found : 0.0;
  return res;
}

}  // namespace radar::experiments
