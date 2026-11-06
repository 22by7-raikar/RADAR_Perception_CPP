#pragma once

// Simulation of real-world impairments the ideal point-target model leaves
// out: mutual interference from another FMCW radar, distributed clutter, and
// multipath reflections off a long flat reflector (guardrail, wall).

#include <cstdint>
#include <vector>

#include "radar/config.hpp"
#include "radar/geometry.hpp"
#include "radar/simulator.hpp"
#include "radar/types.hpp"

namespace radar {

/// Another FMCW radar transmitting towards this one (direct path).
///
/// Its chirps sweep at a different slope, so whenever its instantaneous
/// frequency crosses this radar's ramp the difference frequency sweeps
/// through the receiver's IF band and leaves a short, broadband chirp-like
/// burst in the ADC samples. After the range FFT that energy spreads over all
/// range bins and raises the noise floor.
struct Interferer {
  double start_frequency_hz = 76.6e9;
  double slope_hz_per_s = 3.0e13;  ///< this radar: 1.25e13
  double ramp_duration_s = 30e-6;
  double period_s = 37e-6;
  double time_offset_s = 0.0;  ///< start of its first ramp, frame time
  double inr_db = 30.0;        ///< interference-to-noise ratio per ADC sample
  double azimuth_rad = 0.0;
};

/// Adds the interferer's IF bursts to a simulated frame. The IF band passes
/// beat frequencies in [0, sample_rate), like the target echoes.
void add_interference(const RadarConfig& radar, const Interferer& interferer, DataCube& frame);

/// Point-scatterer clutter patch (vegetation, rough ground): `count` static
/// scatterers uniformly in a range-azimuth sector with internal motion (wind)
/// of standard deviation `velocity_spread_mps`.
struct ClutterPatch {
  double range_min_m = 40.0;
  double range_max_m = 80.0;
  double azimuth_min_rad = -0.9;
  double azimuth_max_rad = -0.4;
  std::size_t count = 400;
  double rcs_dbsm = -12.0;  ///< per scatterer
  double velocity_spread_mps = 0.3;
};

std::vector<Target> clutter_targets(const ClutterPatch& patch, std::uint32_t seed);

/// A long specular reflector parallel to the sensor's x axis at y = `y_m`
/// (sensor frame), such as a guardrail or a tunnel wall.
struct PlanarReflector {
  double y_m = -4.0;
  double loss_db = 3.0;  ///< per bounce
};

/// A reflector in the sensor frame with velocity relative to the sensor.
struct MovingPoint {
  Vec2 position;
  Vec2 velocity;
  double rcs_dbsm = 10.0;
};

/// Direct-path target for a point.
Target direct_target(const MovingPoint& p);

/// Multipath echoes of a point via the reflector (single-transmitter radar):
///  - double bounce (out and back via the reflector): appears at the point's
///    mirror image across the reflector, 2 x loss;
///  - two single-bounce paths (out direct, back via the reflector, and the
///    reverse) with the same delay: range and Doppler halfway between the
///    point and its image. One arrives from the point's direction, one from
///    the image's; 1 x loss each.
/// Returns nothing if the point is on the far side of the reflector.
std::vector<Target> multipath_targets(const MovingPoint& p, const PlanarReflector& wall);

}  // namespace radar
