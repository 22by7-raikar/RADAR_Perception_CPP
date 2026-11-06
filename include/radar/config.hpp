#pragma once

#include <cstddef>
#include <vector>

namespace radar {

inline constexpr double kSpeedOfLight = 299'792'458.0;
inline constexpr double kPi = 3.14159265358979323846;

/// Antenna position in wavelengths: y along the azimuth axis (+y left), z up.
struct ArrayElement {
  double y = 0.0;
  double z = 0.0;
};

/// FMCW front-end parameters and the resolution limits they imply.
///
/// The ADC samples the full chirp ramp, so the sampled bandwidth equals
/// `bandwidth_hz`. Data is complex (IQ), so all `samples_per_chirp` range bins
/// map to non-negative beat frequencies.
///
/// MIMO: with `num_tx > 1` the transmitters take turns (time-division
/// multiplexing). One "loop" is `num_tx` consecutive chirps, one per Tx in
/// the order of `tx_positions`; `chirps_per_frame` counts loops, i.e. chirps
/// per transmitter. Every Tx/Rx pair forms a virtual element at
/// tx_position + rx_position, so frames carry `num_virtual()` channels,
/// ordered channel = tx * num_rx + rx. Receivers form a uniform line along +y.
struct RadarConfig {
  double carrier_hz = 77e9;
  double bandwidth_hz = 150e6;
  double chirp_duration_s = 12e-6;  ///< ADC sampling window per chirp
  double idle_time_s = 2e-6;        ///< gap between chirps
  std::size_t samples_per_chirp = 256;
  std::size_t chirps_per_frame = 128;  ///< chirps per transmitter (Doppler FFT size)
  std::size_t num_rx = 8;
  double rx_spacing_wavelengths = 0.5;  ///< uniform linear array along +y
  double noise_power = 1.0;             ///< complex noise power per ADC sample
  std::size_t num_tx = 1;
  /// Transmitter positions in transmit order. Empty means a single Tx at the
  /// origin (num_tx must then be 1).
  std::vector<ArrayElement> tx_positions;

  double wavelength() const { return kSpeedOfLight / carrier_hz; }
  double slope() const { return bandwidth_hz / chirp_duration_s; }
  double chirp_period() const { return chirp_duration_s + idle_time_s; }
  /// Time between two chirps of the same transmitter.
  double loop_period() const { return static_cast<double>(num_tx) * chirp_period(); }
  double sample_rate() const { return static_cast<double>(samples_per_chirp) / chirp_duration_s; }
  double frame_duration() const { return static_cast<double>(chirps_per_frame) * loop_period(); }

  std::size_t num_virtual() const { return num_tx * num_rx; }
  ArrayElement tx_position(std::size_t tx) const {
    return tx_positions.empty() ? ArrayElement{} : tx_positions[tx];
  }
  ArrayElement rx_position(std::size_t rx) const {
    return {static_cast<double>(rx) * rx_spacing_wavelengths, 0.0};
  }
  /// Position of virtual channel `ch` (= tx * num_rx + rx).
  ArrayElement virtual_position(std::size_t ch) const {
    const ArrayElement t = tx_position(ch / num_rx);
    const ArrayElement r = rx_position(ch % num_rx);
    return {t.y + r.y, t.z + r.z};
  }

  double range_resolution() const { return kSpeedOfLight / (2.0 * bandwidth_hz); }
  double max_range() const { return range_resolution() * static_cast<double>(samples_per_chirp); }
  double velocity_resolution() const {
    return wavelength() / (2.0 * static_cast<double>(chirps_per_frame) * loop_period());
  }
  /// Unambiguous velocity of the Doppler FFT is +/- this value. TDM shrinks it
  /// by num_tx, because each transmitter only chirps every loop_period().
  double max_velocity() const { return wavelength() / (4.0 * loop_period()); }
  /// Velocity span recoverable by TDM phase-hypothesis unfolding (the
  /// single-transmitter limit, set by the chirp period).
  double max_unfolded_velocity() const { return wavelength() / (4.0 * chirp_period()); }

  /// Throws std::invalid_argument if the configuration is not usable.
  void validate() const;

  /// 3 Tx x 4 Rx TDM-MIMO layout common on single-chip 77 GHz radars: Tx0 and
  /// Tx2 form an 8-element azimuth row at lambda/2 spacing; Tx1 is raised by
  /// lambda/2 and adds a 4-element row above the middle of it for elevation.
  static RadarConfig tdm_mimo_3tx4rx();
};

}  // namespace radar
