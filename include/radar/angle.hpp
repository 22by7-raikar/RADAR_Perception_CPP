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
  double peak_power = 0.0;
};

/// Azimuth from the phase progression across a uniform linear array,
/// via a zero-padded FFT over channels with parabolic peak interpolation.
class AngleEstimator {
 public:
  AngleEstimator(const RadarConfig& radar, const AngleConfig& cfg = {});

  /// `spectrum` is the range-Doppler cube ([rx][doppler][range]).
  AngleEstimate estimate(const DataCube& spectrum, std::size_t range_bin, std::size_t doppler_bin);

  /// All resolvable azimuth peaks in the cell, strongest first.
  void estimate_peaks(const DataCube& spectrum, std::size_t range_bin, std::size_t doppler_bin,
                      std::vector<AngleEstimate>& peaks);

 private:
  RadarConfig radar_;
  AngleConfig cfg_;
  std::vector<float> taper_;
  std::vector<cfloat> buffer_;
  std::vector<float> power_;
  std::vector<std::size_t> maxima_;

  double bin_to_azimuth(double bin) const;
  void compute_spectrum(const DataCube& spectrum, std::size_t range_bin, std::size_t doppler_bin);
};

}  // namespace radar
