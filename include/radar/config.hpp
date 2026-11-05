#pragma once

#include <cstddef>

namespace radar {

inline constexpr double kSpeedOfLight = 299'792'458.0;
inline constexpr double kPi = 3.14159265358979323846;

/// FMCW front-end parameters and the resolution limits they imply.
///
/// The ADC samples the full chirp ramp, so the sampled bandwidth equals
/// `bandwidth_hz`. Data is complex (IQ), so all `samples_per_chirp` range bins
/// map to non-negative beat frequencies.
struct RadarConfig {
  double carrier_hz = 77e9;
  double bandwidth_hz = 150e6;
  double chirp_duration_s = 12e-6;  ///< ADC sampling window per chirp
  double idle_time_s = 2e-6;        ///< gap between chirps
  std::size_t samples_per_chirp = 256;
  std::size_t chirps_per_frame = 128;
  std::size_t num_rx = 8;
  double rx_spacing_wavelengths = 0.5;  ///< uniform linear array along +y
  double noise_power = 1.0;             ///< complex noise power per ADC sample

  double wavelength() const { return kSpeedOfLight / carrier_hz; }
  double slope() const { return bandwidth_hz / chirp_duration_s; }
  double chirp_period() const { return chirp_duration_s + idle_time_s; }
  double sample_rate() const { return static_cast<double>(samples_per_chirp) / chirp_duration_s; }
  double frame_duration() const { return static_cast<double>(chirps_per_frame) * chirp_period(); }

  double range_resolution() const { return kSpeedOfLight / (2.0 * bandwidth_hz); }
  double max_range() const { return range_resolution() * static_cast<double>(samples_per_chirp); }
  double velocity_resolution() const {
    return wavelength() / (2.0 * static_cast<double>(chirps_per_frame) * chirp_period());
  }
  /// Unambiguous velocity is +/- this value.
  double max_velocity() const { return wavelength() / (4.0 * chirp_period()); }

  /// Throws std::invalid_argument if the configuration is not usable.
  void validate() const;
};

}  // namespace radar
