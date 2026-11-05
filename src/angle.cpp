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

void AngleEstimator::compute_spectrum(const DataCube& spectrum, std::size_t range_bin,
                                      std::size_t doppler_bin) {
  const std::size_t n = cfg_.fft_size;
  std::fill(buffer_.begin(), buffer_.end(), cfloat{});
  for (std::size_t k = 0; k < radar_.num_rx; ++k) {
    buffer_[k] = spectrum(k, doppler_bin, range_bin) * taper_[k];
  }
  const auto elem = static_cast<std::ptrdiff_t>(sizeof(cfloat));
  pocketfft::c2c({n}, {elem}, {elem}, {0}, pocketfft::FORWARD, buffer_.data(), buffer_.data(),
                 1.0f);
  for (std::size_t i = 0; i < n; ++i) power_[i] = squared_magnitude(buffer_[i]);
}

double AngleEstimator::bin_to_azimuth(double bin) const {
  const auto n = static_cast<double>(cfg_.fft_size);
  if (bin >= n / 2.0) bin -= n;
  // Element k carries phase -2 pi k (d / lambda) sin(az), so the FFT peaks at
  // bin / n = -(d / lambda) sin(az).
  const double sin_az = std::clamp(-bin / (n * radar_.rx_spacing_wavelengths), -1.0, 1.0);
  return std::asin(sin_az);
}

void AngleEstimator::estimate_peaks(const DataCube& spectrum, std::size_t range_bin,
                                    std::size_t doppler_bin, std::vector<AngleEstimate>& peaks) {
  peaks.clear();
  compute_spectrum(spectrum, range_bin, doppler_bin);
  const std::size_t n = cfg_.fft_size;
  maxima_.clear();
  for (std::size_t i = 0; i < n; ++i) {
    const float l = power_[(i + n - 1) % n], c = power_[i], r = power_[(i + 1) % n];
    if (c > l && c >= r) maxima_.push_back(i);
  }
  if (maxima_.empty()) maxima_.push_back(0);  // flat spectrum
  std::sort(maxima_.begin(), maxima_.end(),
            [this](std::size_t a, std::size_t b) { return power_[a] > power_[b]; });

  const double floor =
      static_cast<double>(power_[maxima_.front()]) * std::pow(10.0, cfg_.secondary_peak_db / 10.0);
  for (std::size_t i : maxima_) {
    if (peaks.size() >= std::max<std::size_t>(1, cfg_.max_peaks)) break;
    if (!peaks.empty() && static_cast<double>(power_[i]) < floor) break;
    const double delta = parabolic_offset(power_[(i + n - 1) % n], power_[i], power_[(i + 1) % n]);
    peaks.push_back(
        {bin_to_azimuth(static_cast<double>(i) + delta), static_cast<double>(power_[i])});
  }
}

AngleEstimate AngleEstimator::estimate(const DataCube& spectrum, std::size_t range_bin,
                                       std::size_t doppler_bin) {
  std::vector<AngleEstimate> peaks;
  estimate_peaks(spectrum, range_bin, doppler_bin, peaks);
  return peaks.front();
}

}  // namespace radar
