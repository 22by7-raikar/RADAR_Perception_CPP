#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "radar/config.hpp"
#include "radar/types.hpp"

namespace radar {

/// Point target with constant radial velocity.
struct Target {
  double range_m = 0.0;
  double velocity_mps = 0.0;  ///< range rate, negative = approaching
  double azimuth_rad = 0.0;   ///< counter-clockwise from boresight (+x)
  double rcs_dbsm = 10.0;
};

struct SimulatorConfig {
  /// Per-sample SNR of a reference target; other targets scale with the
  /// radar equation (SNR ~ RCS / R^4).
  double reference_snr_db = -10.0;
  double reference_range_m = 100.0;
  double reference_rcs_dbsm = 10.0;
  bool add_noise = true;
  std::uint32_t seed = 42;
};

/// Synthesises the dechirped IQ beat signal of an FMCW radar.
///
/// Each sample uses the exact round-trip delay
///   tau = 2 (R0 + v t) / c - k d sin(az) / c   (t from the frame start)
/// and the dechirped phase 2 pi (fc tau + S t_fast tau - S tau^2 / 2),
/// so range migration and range-Doppler coupling are included.
class Simulator {
 public:
  Simulator(const RadarConfig& radar, const SimulatorConfig& sim = {});

  /// Simulates one frame. Target states are taken at the start of the first
  /// chirp; targets keep moving at constant velocity through the frame.
  void simulate_frame(const std::vector<Target>& targets, DataCube& out);
  DataCube simulate_frame(const std::vector<Target>& targets);

  /// Linear amplitude of a target's echo per ADC sample.
  double echo_amplitude(const Target& target) const;

  const RadarConfig& radar() const { return radar_; }

 private:
  RadarConfig radar_;
  SimulatorConfig sim_;
  std::mt19937 rng_;
};

}  // namespace radar
