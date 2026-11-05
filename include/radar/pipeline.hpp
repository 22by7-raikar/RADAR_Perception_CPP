#pragma once

#include <vector>

#include "radar/angle.hpp"
#include "radar/cfar.hpp"
#include "radar/config.hpp"
#include "radar/range_doppler.hpp"
#include "radar/types.hpp"

namespace radar {

struct PipelineConfig {
  RangeDopplerConfig range_doppler;
  CfarConfig cfar;  ///< num_looks is overridden with the number of channels
  AngleConfig angle;
  bool peak_grouping = true;
  double min_range_m = 1.0;  ///< drop detections closer than this (DC / leakage)
};

struct StageTimings {
  double range_doppler_ms = 0.0;
  double integration_ms = 0.0;
  double cfar_ms = 0.0;
  double peak_grouping_ms = 0.0;
  double angle_ms = 0.0;
  double total_ms = 0.0;
};

struct FrameResult {
  std::vector<Detection> detections;
  std::size_t cfar_cells = 0;  ///< cells above threshold before peak grouping
  StageTimings timings;
};

/// ADC frame -> point cloud.
///
/// All working buffers are owned by the pipeline and reused, so steady-state
/// processing does not allocate (beyond output growth).
class RadarPipeline {
 public:
  RadarPipeline(const RadarConfig& radar, const PipelineConfig& cfg = {});

  void process(const DataCube& adc, FrameResult& result);

  const Grid2D<float>& power_map() const { return power_; }
  const DataCube& spectrum() const { return spectrum_; }
  const RadarConfig& radar() const { return radar_; }
  const RangeDopplerProcessor& range_doppler() const { return rd_; }

 private:
  RadarConfig radar_;
  PipelineConfig cfg_;
  RangeDopplerProcessor rd_;
  CaCfar2D cfar_;
  AngleEstimator angle_;

  DataCube spectrum_;
  Grid2D<float> power_;
  std::vector<CfarCell> cells_;
  std::vector<CfarCell> peaks_;
};

}  // namespace radar
