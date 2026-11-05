#include "radar/angle.hpp"

#include <pocketfft_hdronly.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "radar/interpolation.hpp"

namespace radar {

AngleEstimator::AngleEstimator(const RadarConfig& radar, const AngleConfig& cfg)
    : radar_(radar), cfg_(cfg) {
  radar_.validate();
  if (cfg_.fft_size < radar_.num_rx) {
    throw std::invalid_argument("AngleEstimator: fft_size must be >= num_rx");
  }
  taper_ = make_window(cfg_.window, radar_.num_rx);
  buffer_.resize(cfg_.fft_size);
  power_.resize(cfg_.fft_size);
}

AngleEstimate AngleEstimator::estimate(const DataCube& spectrum, std::size_t range_bin,
                                       std::size_t doppler_bin) {
  const std::size_t n = cfg_.fft_size;
  std::fill(buffer_.begin(), buffer_.end(), cfloat{});
  for (std::size_t k = 0; k < radar_.num_rx; ++k) {
    buffer_[k] = spectrum(k, doppler_bin, range_bin) * taper_[k];
  }
  const auto elem = static_cast<std::ptrdiff_t>(sizeof(cfloat));
  pocketfft::c2c({n}, {elem}, {elem}, {0}, pocketfft::FORWARD, buffer_.data(), buffer_.data(),
                 1.0f);

  std::size_t peak = 0;
  for (std::size_t i = 0; i < n; ++i) {
    power_[i] = std::norm(buffer_[i]);
    if (power_[i] > power_[peak]) peak = i;
  }
  const double delta =
      parabolic_offset(power_[(peak + n - 1) % n], power_[peak], power_[(peak + 1) % n]);
  double bin = static_cast<double>(peak) + delta;
  if (bin >= static_cast<double>(n) / 2.0) bin -= static_cast<double>(n);

  // Element k carries phase -2 pi k (d / lambda) sin(az), so the FFT peaks at
  // bin / n = -(d / lambda) sin(az).
  const double sin_az =
      std::clamp(-bin / (static_cast<double>(n) * radar_.rx_spacing_wavelengths), -1.0, 1.0);
  return {std::asin(sin_az), static_cast<double>(power_[peak])};
}

}  // namespace radar
