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

 private:
  RadarConfig radar_;
  AngleConfig cfg_;
  std::vector<float> taper_;
  std::vector<cfloat> buffer_;
  std::vector<float> power_;
};

}  // namespace radar
