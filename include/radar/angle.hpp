#pragma once

#include <cstddef>
#include <vector>

#include "radar/config.hpp"
#include "radar/types.hpp"
#include "radar/window.hpp"

namespace radar {

struct AngleConfig {
  std::size_t fft_size = 64;  ///< zero-padded angle FFT length
  WindowType window = WindowType::kRectangular;
  /// Up to this many azimuth peaks are reported per range-Doppler cell, so two
  /// reflectors that share a cell (e.g. a slow pedestrian next to a static
  /// post) are both detected when the array can separate them.
  std::size_t max_peaks = 2;
  /// Secondary peaks must be within this many dB of the strongest one. Kept
  /// above the -13 dB first sidelobe of the rectangular-tapered array.
  double secondary_peak_db = -9.0;
};

struct AngleEstimate {
  double azimuth_rad = 0.0;
  double elevation_rad = 0.0;  ///< 0 unless the array has an elevation row
  double peak_power = 0.0;
};

/// The parts of a virtual array the angle estimator uses.
///
/// The azimuth row is every element at z = 0, sorted by y; it must be a
/// uniform line. Elevation pairs are elements above the row (all at the same
/// height dz) matched with the row element at the same y.
struct VirtualArray {
  std::vector<std::size_t> azimuth_row;  ///< channel indices, increasing y
  double azimuth_spacing = 0.5;          ///< wavelengths
  struct ElevationPair {
    std::size_t row;     ///< channel in the azimuth row
    std::size_t raised;  ///< channel dz above it
    double y;            ///< shared y position, wavelengths
  };
  std::vector<ElevationPair> elevation_pairs;
  double elevation_spacing = 0.0;  ///< dz in wavelengths
  bool has_elevation() const { return !elevation_pairs.empty(); }
};

/// Throws std::invalid_argument if the azimuth row is not a uniform line.
VirtualArray analyse_virtual_array(const RadarConfig& radar);

/// Azimuth from the phase progression across the azimuth row (zero-padded FFT
/// with parabolic peak interpolation); elevation, when the array has raised
/// elements, from the phase between the raised row and the elements below it
/// after beamforming both at the azimuth peak.
///
/// Element (y, z) of a target at direction u = cos(el) sin(az), w = sin(el)
/// carries phase -2 pi (y u + z w). The FFT measures u; elevation then
/// corrects the azimuth to asin(u / cos(el)).
class AngleEstimator {
 public:
  AngleEstimator(const RadarConfig& radar, const AngleConfig& cfg = {});

  /// `spectrum` is the range-Doppler cube ([channel][doppler][range]).
  AngleEstimate estimate(const DataCube& spectrum, std::size_t range_bin, std::size_t doppler_bin);

  /// All resolvable azimuth peaks in the cell, strongest first.
  void estimate_peaks(const DataCube& spectrum, std::size_t range_bin, std::size_t doppler_bin,
                      std::vector<AngleEstimate>& peaks);

  /// Same, from one snapshot of `num_virtual()` channels (after any TDM
  /// motion compensation).
  void estimate_peaks(const cfloat* snapshot, std::vector<AngleEstimate>& peaks);

  const VirtualArray& array() const { return array_; }

 private:
  RadarConfig radar_;
  AngleConfig cfg_;
  VirtualArray array_;
  std::vector<float> taper_;
  std::vector<cfloat> snapshot_;
  std::vector<cfloat> buffer_;
  std::vector<float> power_;
  std::vector<std::size_t> maxima_;

  double bin_to_u(double bin) const;
  double estimate_elevation(const cfloat* snapshot, double u) const;
};

}  // namespace radar
