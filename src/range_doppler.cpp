#include "radar/range_doppler.hpp"

#include <pocketfft_hdronly.h>

#include <algorithm>
#include <stdexcept>

namespace radar {

RangeDopplerProcessor::RangeDopplerProcessor(const RadarConfig& radar,
                                             const RangeDopplerConfig& cfg)
    : radar_(radar) {
  radar_.validate();
  const std::size_t nd = radar_.chirps_per_frame;
  const std::size_t nr = radar_.samples_per_chirp;
  const std::vector<float> wr = make_window(cfg.range_window, nr);
  const std::vector<float> wd = make_window(cfg.doppler_window, nd);

  taper_.resize(nd * nr);
  for (std::size_t m = 0; m < nd; ++m) {
    const float sign = (m % 2 == 0) ? 1.0f : -1.0f;  // shifts zero Doppler to bin nd/2
    for (std::size_t n = 0; n < nr; ++n) taper_[m * nr + n] = sign * wd[m] * wr[n];
  }
}

void RangeDopplerProcessor::process(const DataCube& adc, DataCube& spectrum) const {
  const std::size_t nrx = radar_.num_rx;
  const std::size_t nd = radar_.chirps_per_frame;
  const std::size_t nr = radar_.samples_per_chirp;
  if (adc.num_rx() != nrx || adc.num_chirps() != nd || adc.num_samples() != nr) {
    throw std::invalid_argument("RangeDopplerProcessor: frame shape does not match config");
  }
  if (spectrum.num_rx() != nrx || spectrum.num_chirps() != nd || spectrum.num_samples() != nr) {
    spectrum.resize(nrx, nd, nr);
  }

  const std::size_t plane = nd * nr;
  for (std::size_t k = 0; k < nrx; ++k) {
    const cfloat* src = adc.data() + k * plane;
    cfloat* dst = spectrum.data() + k * plane;
    for (std::size_t i = 0; i < plane; ++i) dst[i] = src[i] * taper_[i];
  }

  const pocketfft::shape_t shape{nrx, nd, nr};
  const auto elem = static_cast<std::ptrdiff_t>(sizeof(cfloat));
  const pocketfft::stride_t stride{elem * static_cast<std::ptrdiff_t>(plane),
                                   elem * static_cast<std::ptrdiff_t>(nr), elem};
  pocketfft::c2c(shape, stride, stride, {1, 2}, pocketfft::FORWARD, spectrum.data(),
                 spectrum.data(), 1.0f);
}

double RangeDopplerProcessor::bin_to_velocity(double doppler_bin) const {
  return (doppler_bin - static_cast<double>(radar_.chirps_per_frame / 2)) *
         radar_.velocity_resolution();
}

double RangeDopplerProcessor::bin_to_range(double range_bin) const {
  return range_bin * radar_.range_resolution();
}

void integrate_power(const DataCube& spectrum, Grid2D<float>& power) {
  const std::size_t nrx = spectrum.num_rx();
  const std::size_t nd = spectrum.num_chirps();
  const std::size_t nr = spectrum.num_samples();
  if (power.rows() != nr || power.cols() != nd) power.resize(nr, nd);

  // Accumulate each Doppler row contiguously across channels, then write the
  // column once. Avoids std::norm, which libstdc++ implements via hypot().
  thread_local std::vector<float> acc;
  acc.resize(nr);
  for (std::size_t d = 0; d < nd; ++d) {
    std::fill(acc.begin(), acc.end(), 0.0f);
    for (std::size_t k = 0; k < nrx; ++k) {
      const cfloat* row = &spectrum(k, d, 0);
      for (std::size_t r = 0; r < nr; ++r) acc[r] += squared_magnitude(row[r]);
    }
    for (std::size_t r = 0; r < nr; ++r) power(r, d) = acc[r];
  }
}

}  // namespace radar
