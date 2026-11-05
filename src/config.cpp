#include "radar/config.hpp"

#include <stdexcept>
#include <string>

namespace radar {

namespace {
void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("RadarConfig: ") + what);
}
}  // namespace

void RadarConfig::validate() const {
  require(carrier_hz > 0.0, "carrier_hz must be positive");
  require(bandwidth_hz > 0.0, "bandwidth_hz must be positive");
  require(chirp_duration_s > 0.0, "chirp_duration_s must be positive");
  require(idle_time_s >= 0.0, "idle_time_s must be non-negative");
  require(samples_per_chirp >= 8, "samples_per_chirp must be >= 8");
  require(chirps_per_frame >= 8 && chirps_per_frame % 2 == 0,
          "chirps_per_frame must be even and >= 8");
  require(num_rx >= 1, "num_rx must be >= 1");
  require(rx_spacing_wavelengths > 0.0, "rx_spacing_wavelengths must be positive");
  require(noise_power >= 0.0, "noise_power must be non-negative");
}

}  // namespace radar
