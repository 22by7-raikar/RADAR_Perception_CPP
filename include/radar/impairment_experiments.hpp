#pragma once

// Fixed-seed experiments for the impairments and their mitigations, shared by
// impairments_demo (figures, report) and the regression suite (bounds).

#include <cstdint>
#include <vector>

#include "radar/mitigation.hpp"
#include "radar/simulator.hpp"
#include "radar/types.hpp"

namespace radar::experiments {

struct DetectionStats {
  double detection_rate = 0.0;
  double false_alarms_per_frame = 0.0;
};

// --- 1. Interference from another FMCW radar ---------------------------------

struct InterferenceResult {
  DetectionStats clean, interfered, mitigated;
  double floor_rise_db = 0.0;           ///< interfered vs. clean noise floor
  double residual_floor_db = 0.0;       ///< after blanking vs. clean
  double blanked_fraction = 0.0;        ///< of all ADC (chirp, sample) slots
  double weak_target_rate_clean = 0.0;  ///< targets <= 0 dBsm beyond 100 m
  double weak_target_rate_interfered = 0.0;
  double weak_target_rate_mitigated = 0.0;
  double frame_ms_plain = 0.0;      ///< full pipeline on the interfered frame
  double frame_ms_mitigated = 0.0;  ///< same, with detection + blanking + IMAT

  // Last frame, for figures.
  std::vector<float> chirp_envelope;  ///< summed over Rx, one chirp with a burst
  std::vector<float> chirp_residual;  ///< Rx 0 with echoes notched out (detector input)
  std::vector<char> chirp_blanked;
  std::vector<float> range_profile_clean_db;  ///< mean over Doppler
  std::vector<float> range_profile_interfered_db;
  std::vector<float> range_profile_mitigated_db;
};

InterferenceResult run_interference(int frames, std::uint32_t seed);

// --- 2. Clutter and target masking: CA- vs OS-CFAR ---------------------------

struct CfarComparison {
  double clutter_edge_rate = 0.0;  ///< slow pedestrian 3 m in front of vegetation
  double masked_rate = 0.0;        ///< weak target 4 m behind a strong one
  double reference_rate = 0.0;     ///< isolated target of the same strength
  double clutter_detections_per_frame = 0.0;
  double false_alarms_per_frame = 0.0;  ///< outside targets and clutter
  double cfar_ms = 0.0;
};

struct ClutterResult {
  CfarComparison ca, os;

  // Last frame, for figures: power map (range x Doppler, dB) and detections.
  std::size_t rows = 0, cols = 0;
  std::vector<float> power_db;
  std::vector<Detection> ca_detections, os_detections;
  std::vector<Target> truth;  ///< pedestrian, car, motorbike, reference
};

ClutterResult run_clutter(int frames, std::uint32_t seed);

// --- 3. Multipath ghosts off a guardrail -------------------------------------

enum class Source { kReal, kRail, kGhostMirror, kGhostDirect };

struct MultipathResult {
  double mirror_ghosts_per_frame = 0.0;       ///< double bounce + image-side single bounce
  double mirror_ghosts_left_per_frame = 0.0;  ///< after the filter
  double direct_ghosts_per_frame = 0.0;       ///< single bounce arriving along the real bearing
  double real_flagged_per_frame = 0.0;        ///< real moving detections wrongly removed
  double real_moving_per_frame = 0.0;
  double rail_offset_error_m = 0.0;  ///< mean |fitted - true| at x = 0
  double rail_found_fraction = 0.0;

  // Last frame, for figures.
  std::vector<Detection> detections;
  std::vector<Source> sources;
  std::vector<char> flagged;
  ReflectorLine rail;
};

MultipathResult run_multipath(int frames, std::uint32_t seed);

}  // namespace radar::experiments
