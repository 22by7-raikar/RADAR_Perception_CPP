#include "radar/io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "radar/config.hpp"

namespace radar {

namespace {

std::ofstream open_out(const std::string& path, bool append) {
  std::ofstream f(path, append ? std::ios::app : std::ios::trunc);
  if (!f) throw std::runtime_error("Cannot open " + path + " for writing");
  return f;
}

}  // namespace

std::vector<Target> load_targets_csv(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("Cannot open scenario file " + path);

  std::vector<Target> targets;
  std::string line;
  std::size_t line_no = 0;
  while (std::getline(f, line)) {
    ++line_no;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    if (std::isalpha(static_cast<unsigned char>(line[first]))) continue;  // header

    std::stringstream ss(line);
    std::string field;
    std::vector<double> values;
    while (std::getline(ss, field, ',')) {
      try {
        values.push_back(std::stod(field));
      } catch (const std::exception&) {
        throw std::runtime_error(path + ":" + std::to_string(line_no) + ": bad number '" + field +
                                 "'");
      }
    }
    if (values.size() < 3) {
      throw std::runtime_error(path + ":" + std::to_string(line_no) + ": expected >= 3 columns");
    }
    Target t;
    t.range_m = values[0];
    t.velocity_mps = values[1];
    t.azimuth_rad = values[2] * kPi / 180.0;
    if (values.size() > 3) t.rcs_dbsm = values[3];
    targets.push_back(t);
  }
  return targets;
}

void write_power_map_csv(const std::string& path, const Grid2D<float>& power) {
  auto f = open_out(path, false);
  for (std::size_t r = 0; r < power.rows(); ++r) {
    for (std::size_t d = 0; d < power.cols(); ++d) {
      if (d) f << ',';
      f << 10.0 * std::log10(std::max(1e-30, static_cast<double>(power(r, d))));
    }
    f << '\n';
  }
}

void write_detections_csv(const std::string& path, const std::vector<Detection>& detections,
                          std::size_t frame_index, bool append) {
  auto f = open_out(path, append);
  if (!append) f << "frame,range_m,velocity_mps,azimuth_deg,x_m,y_m,snr_db,range_bin,doppler_bin\n";
  for (const Detection& d : detections) {
    f << frame_index << ',' << d.range_m << ',' << d.velocity_mps << ','
      << d.azimuth_rad * 180.0 / kPi << ',' << d.x_m << ',' << d.y_m << ',' << d.snr_db << ','
      << d.range_bin << ',' << d.doppler_bin << '\n';
  }
}

void write_targets_csv(const std::string& path, const std::vector<Target>& targets,
                       std::size_t frame_index, bool append) {
  auto f = open_out(path, append);
  if (!append) f << "frame,range_m,velocity_mps,azimuth_deg,rcs_dbsm\n";
  for (const Target& t : targets) {
    f << frame_index << ',' << t.range_m << ',' << t.velocity_mps << ','
      << t.azimuth_rad * 180.0 / kPi << ',' << t.rcs_dbsm << '\n';
  }
}

}  // namespace radar
