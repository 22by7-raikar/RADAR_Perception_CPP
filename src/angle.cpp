#include "radar/angle.hpp"

#include <pocketfft_hdronly.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "radar/interpolation.hpp"

namespace radar {

VirtualArray analyse_virtual_array(const RadarConfig& radar) {
  radar.validate();
  constexpr double kTol = 1e-6;
  const std::size_t n = radar.num_virtual();
  VirtualArray va;

  // Azimuth row: z = 0, sorted by y, duplicates (overlapping elements) dropped.
  std::vector<std::size_t> row;
  for (std::size_t ch = 0; ch < n; ++ch) {
    if (std::abs(radar.virtual_position(ch).z) < kTol) row.push_back(ch);
  }
  std::stable_sort(row.begin(), row.end(), [&](std::size_t a, std::size_t b) {
    return radar.virtual_position(a).y < radar.virtual_position(b).y;
  });
  for (std::size_t ch : row) {
    if (!va.azimuth_row.empty() && std::abs(radar.virtual_position(ch).y -
                                            radar.virtual_position(va.azimuth_row.back()).y) < kTol)
      continue;
    va.azimuth_row.push_back(ch);
  }
  if (va.azimuth_row.empty()) throw std::invalid_argument("Virtual array has no element at z = 0");
  if (va.azimuth_row.size() > 1) {
    const double y0 = radar.virtual_position(va.azimuth_row[0]).y;
    va.azimuth_spacing = radar.virtual_position(va.azimuth_row[1]).y - y0;
    for (std::size_t i = 1; i < va.azimuth_row.size(); ++i) {
      const double step = radar.virtual_position(va.azimuth_row[i]).y -
                          radar.virtual_position(va.azimuth_row[i - 1]).y;
      if (std::abs(step - va.azimuth_spacing) > kTol) {
        throw std::invalid_argument("Virtual azimuth row is not uniformly spaced");
      }
    }
  } else {
    va.azimuth_spacing = radar.rx_spacing_wavelengths;
  }

  // Elevation pairs: raised elements with a row element at the same y.
  for (std::size_t ch = 0; ch < n; ++ch) {
    const ArrayElement p = radar.virtual_position(ch);
    if (std::abs(p.z) < kTol) continue;
    if (va.elevation_spacing == 0.0) va.elevation_spacing = p.z;
    if (std::abs(p.z - va.elevation_spacing) > kTol) {
      throw std::invalid_argument("Raised virtual elements must share one height");
    }
    for (std::size_t r : va.azimuth_row) {
      if (std::abs(radar.virtual_position(r).y - p.y) < kTol) {
        va.elevation_pairs.push_back({r, ch, p.y});
        break;
      }
    }
  }
  return va;
}

AngleEstimator::AngleEstimator(const RadarConfig& radar, const AngleConfig& cfg)
    : radar_(radar), cfg_(cfg), array_(analyse_virtual_array(radar)) {
  if (cfg_.fft_size < array_.azimuth_row.size()) {
    throw std::invalid_argument("AngleEstimator: fft_size must be >= the azimuth row length");
  }
  taper_ = make_window(cfg_.window, array_.azimuth_row.size());
  snapshot_.resize(radar_.num_virtual());
  buffer_.resize(cfg_.fft_size);
  power_.resize(cfg_.fft_size);
}

double AngleEstimator::bin_to_u(double bin) const {
  const auto n = static_cast<double>(cfg_.fft_size);
  if (bin >= n / 2.0) bin -= n;
  // Row element i carries phase -2 pi i d u, so the FFT peaks at bin / n = -d u.
  return std::clamp(-bin / (n * array_.azimuth_spacing), -1.0, 1.0);
}

double AngleEstimator::estimate_elevation(const cfloat* snapshot, double u) const {
  // Beamform the raised elements and the row elements below them at u; the
  // phase between the two beams is -2 pi dz w.
  std::complex<double> low{}, high{};
  for (const auto& p : array_.elevation_pairs) {
    const std::complex<double> steer = std::polar(1.0, 2.0 * kPi * p.y * u);
    low += std::complex<double>(snapshot[p.row]) * steer;
    high += std::complex<double>(snapshot[p.raised]) * steer;
  }
  const double phase = std::arg(high * std::conj(low));
  const double w = std::clamp(-phase / (2.0 * kPi * array_.elevation_spacing), -1.0, 1.0);
  return std::asin(w);
}

void AngleEstimator::estimate_peaks(const DataCube& spectrum, std::size_t range_bin,
                                    std::size_t doppler_bin, std::vector<AngleEstimate>& peaks) {
  for (std::size_t ch = 0; ch < snapshot_.size(); ++ch) {
    snapshot_[ch] = spectrum(ch, doppler_bin, range_bin);
  }
  estimate_peaks(snapshot_.data(), peaks);
}

void AngleEstimator::estimate_peaks(const cfloat* snapshot, std::vector<AngleEstimate>& peaks) {
  peaks.clear();
  const std::size_t n = cfg_.fft_size;
  std::fill(buffer_.begin(), buffer_.end(), cfloat{});
  for (std::size_t i = 0; i < array_.azimuth_row.size(); ++i) {
    buffer_[i] = snapshot[array_.azimuth_row[i]] * taper_[i];
  }
  const auto elem = static_cast<std::ptrdiff_t>(sizeof(cfloat));
  pocketfft::c2c({n}, {elem}, {elem}, {0}, pocketfft::FORWARD, buffer_.data(), buffer_.data(),
                 1.0f);
  for (std::size_t i = 0; i < n; ++i) power_[i] = squared_magnitude(buffer_[i]);

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
    const double u = bin_to_u(static_cast<double>(i) + delta);
    AngleEstimate a;
    a.peak_power = static_cast<double>(power_[i]);
    if (array_.has_elevation()) {
      a.elevation_rad = estimate_elevation(snapshot, u);
      a.azimuth_rad = std::asin(std::clamp(u / std::cos(a.elevation_rad), -1.0, 1.0));
    } else {
      a.azimuth_rad = std::asin(u);
    }
    peaks.push_back(a);
  }
}

AngleEstimate AngleEstimator::estimate(const DataCube& spectrum, std::size_t range_bin,
                                       std::size_t doppler_bin) {
  std::vector<AngleEstimate> peaks;
  estimate_peaks(spectrum, range_bin, doppler_bin, peaks);
  return peaks.front();
}

}  // namespace radar
