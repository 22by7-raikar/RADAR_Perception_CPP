#pragma once

#include <string>
#include <vector>

#include "radar/simulator.hpp"
#include "radar/types.hpp"

namespace radar {

/// Reads targets from CSV: range_m, velocity_mps, azimuth_deg[, rcs_dbsm].
/// Blank lines, '#' comments and a header row are skipped.
std::vector<Target> load_targets_csv(const std::string& path);

/// Writes the power map in dB with range rows and Doppler columns.
void write_power_map_csv(const std::string& path, const Grid2D<float>& power);

void write_detections_csv(const std::string& path, const std::vector<Detection>& detections,
                          std::size_t frame_index, bool append);

void write_targets_csv(const std::string& path, const std::vector<Target>& targets,
                       std::size_t frame_index, bool append);

}  // namespace radar
