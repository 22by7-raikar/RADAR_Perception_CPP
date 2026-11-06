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
  double elevation_rad = 0.0;  ///< above the sensor's horizontal plane
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
///   tau = 2 (R0 + v t) / c - lambda (y u + z w) / c   (t from the frame start)
/// where (y, z) is the virtual element position in wavelengths,
/// u = cos(el) sin(az) and w = sin(el), and the dechirped phase
/// 2 pi (fc tau + S t_fast tau - S tau^2 / 2). Range migration, range-Doppler
/// coupling and the TDM-MIMO chirp timing (Tx t transmits at offset t * Tc in
/// every loop) are therefore all included.
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
