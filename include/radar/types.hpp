#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

namespace radar {

using cfloat = std::complex<float>;

/// |z|^2 without the hypot() round trip that std::norm takes in libstdc++.
inline float squared_magnitude(const cfloat& z) {
  return z.real() * z.real() + z.imag() * z.imag();
}

/// Dense, row-major 2D array.
template <typename T>
class Grid2D {
 public:
  Grid2D() = default;
  Grid2D(std::size_t rows, std::size_t cols, T value = T{})
      : rows_(rows), cols_(cols), data_(rows * cols, value) {}

  void resize(std::size_t rows, std::size_t cols) {
    rows_ = rows;
    cols_ = cols;
    data_.resize(rows * cols);
  }
  void fill(T value) { std::fill(data_.begin(), data_.end(), value); }

  T& operator()(std::size_t r, std::size_t c) { return data_[r * cols_ + c]; }
  const T& operator()(std::size_t r, std::size_t c) const { return data_[r * cols_ + c]; }

  std::size_t rows() const { return rows_; }
  std::size_t cols() const { return cols_; }
  std::size_t size() const { return data_.size(); }
  T* data() { return data_.data(); }
  const T* data() const { return data_.data(); }

 private:
  std::size_t rows_ = 0;
  std::size_t cols_ = 0;
  std::vector<T> data_;
};

/// Complex radar data cube, laid out as [rx][chirp][sample].
///
/// The same container holds the raw ADC frame (fast time along the last axis)
/// and the range-Doppler spectrum ([rx][doppler][range]) after processing.
class DataCube {
 public:
  DataCube() = default;
  DataCube(std::size_t num_rx, std::size_t num_chirps, std::size_t num_samples)
      : num_rx_(num_rx),
        num_chirps_(num_chirps),
        num_samples_(num_samples),
        data_(num_rx * num_chirps * num_samples) {}

  void resize(std::size_t num_rx, std::size_t num_chirps, std::size_t num_samples) {
    num_rx_ = num_rx;
    num_chirps_ = num_chirps;
    num_samples_ = num_samples;
    data_.assign(num_rx * num_chirps * num_samples, cfloat{});
  }

  cfloat& operator()(std::size_t rx, std::size_t chirp, std::size_t sample) {
    return data_[(rx * num_chirps_ + chirp) * num_samples_ + sample];
  }
  const cfloat& operator()(std::size_t rx, std::size_t chirp, std::size_t sample) const {
    return data_[(rx * num_chirps_ + chirp) * num_samples_ + sample];
  }

  std::size_t num_rx() const { return num_rx_; }
  std::size_t num_chirps() const { return num_chirps_; }
  std::size_t num_samples() const { return num_samples_; }
  std::size_t size() const { return data_.size(); }
  cfloat* data() { return data_.data(); }
  const cfloat* data() const { return data_.data(); }

 private:
  std::size_t num_rx_ = 0;
  std::size_t num_chirps_ = 0;
  std::size_t num_samples_ = 0;
  std::vector<cfloat> data_;
};

/// A single detected reflection in sensor coordinates.
///
/// Frame convention (ROS REP-103): x forward, y left, azimuth counter-clockwise
/// from +x. Velocity is the radial range rate: negative means approaching.
struct Detection {
  float range_m = 0.0f;
  float velocity_mps = 0.0f;
  float azimuth_rad = 0.0f;
  float x_m = 0.0f;
  float y_m = 0.0f;
  float snr_db = 0.0f;
  std::size_t range_bin = 0;
  std::size_t doppler_bin = 0;
};

}  // namespace radar
