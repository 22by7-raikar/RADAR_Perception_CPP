#include "radar/simulator.hpp"

#include <cmath>
#include <stdexcept>

namespace radar {

Simulator::Simulator(const RadarConfig& radar, const SimulatorConfig& sim)
    : radar_(radar), sim_(sim), rng_(sim.seed) {
  radar_.validate();
}

double Simulator::echo_amplitude(const Target& target) const {
  if (target.range_m <= 0.0) throw std::invalid_argument("Target range must be positive");
  const double range_ratio = sim_.reference_range_m / target.range_m;
  const double snr_db = sim_.reference_snr_db + (target.rcs_dbsm - sim_.reference_rcs_dbsm) +
                        40.0 * std::log10(range_ratio);
  const double noise = radar_.noise_power > 0.0 ? radar_.noise_power : 1.0;
  return std::sqrt(noise * std::pow(10.0, snr_db / 10.0));
}

DataCube Simulator::simulate_frame(const std::vector<Target>& targets) {
  DataCube cube;
  simulate_frame(targets, cube);
  return cube;
}

void Simulator::simulate_frame(const std::vector<Target>& targets, DataCube& out) {
  const std::size_t nrx = radar_.num_rx;
  const std::size_t nd = radar_.chirps_per_frame;
  const std::size_t nr = radar_.samples_per_chirp;
  out.resize(nrx, nd, nr);

  const double c = kSpeedOfLight;
  const double fc = radar_.carrier_hz;
  const double slope = radar_.slope();
  const double dt = 1.0 / radar_.sample_rate();
  const double tp = radar_.chirp_period();
  const double spacing_m = radar_.rx_spacing_wavelengths * radar_.wavelength();

  for (const Target& tgt : targets) {
    const double amp = echo_amplitude(tgt);
    const double sin_az = std::sin(tgt.azimuth_rad);
    for (std::size_t k = 0; k < nrx; ++k) {
      // Element k sits at y = k * d and receives earlier for targets on the left.
      const double array_delay = static_cast<double>(k) * spacing_m * sin_az / c;
      for (std::size_t m = 0; m < nd; ++m) {
        const double t_chirp = static_cast<double>(m) * tp;
        cfloat* row = &out(k, m, 0);
        for (std::size_t n = 0; n < nr; ++n) {
          const double t_fast = static_cast<double>(n) * dt;
          const double tau =
              2.0 * (tgt.range_m + tgt.velocity_mps * (t_chirp + t_fast)) / c - array_delay;
          double cycles = fc * tau + slope * t_fast * tau - 0.5 * slope * tau * tau;
          cycles -= std::floor(cycles);
          const double phase = 2.0 * kPi * cycles;
          row[n] += cfloat(static_cast<float>(amp * std::cos(phase)),
                           static_cast<float>(amp * std::sin(phase)));
        }
      }
    }
  }

  if (sim_.add_noise && radar_.noise_power > 0.0) {
    std::normal_distribution<float> noise(0.0f,
                                          static_cast<float>(std::sqrt(radar_.noise_power / 2.0)));
    cfloat* p = out.data();
    for (std::size_t i = 0; i < out.size(); ++i) p[i] += cfloat(noise(rng_), noise(rng_));
  }
}

}  // namespace radar
