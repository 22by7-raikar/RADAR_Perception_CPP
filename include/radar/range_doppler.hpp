#pragma once

#include <vector>

#include "radar/config.hpp"
#include "radar/types.hpp"
#include "radar/window.hpp"

namespace radar {

struct RangeDopplerConfig {
  WindowType range_window = WindowType::kHann;
  WindowType doppler_window = WindowType::kHann;
};

/// Windowed 2D FFT (range x Doppler) per receive channel.
///
/// Output layout is [rx][doppler][range]. The Doppler axis is centred: bin
/// `chirps_per_frame / 2` is zero velocity. Centring is folded into the window
/// as a (-1)^m modulation, so no separate fftshift pass is needed.
class RangeDopplerProcessor {
 public:
  explicit RangeDopplerProcessor(const RadarConfig& radar, const RangeDopplerConfig& cfg = {});

  void process(const DataCube& adc, DataCube& spectrum) const;

  /// Velocity of a (possibly fractional) Doppler bin.
  double bin_to_velocity(double doppler_bin) const;
  /// Range of a (possibly fractional) range bin.
  double bin_to_range(double range_bin) const;

 private:
  RadarConfig radar_;
  std::vector<float> taper_;  ///< combined [chirp][sample] window incl. Doppler centring
};

/// Non-coherent integration across channels: power(range, doppler) = sum_k |X_k|^2.
void integrate_power(const DataCube& spectrum, Grid2D<float>& power);

}  // namespace radar
