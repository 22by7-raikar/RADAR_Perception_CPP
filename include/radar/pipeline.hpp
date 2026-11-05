#pragma once

#include <cmath>
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
  /// Candidates in neighbouring cells closer than this in azimuth are grouped.
  double grouping_angle_rad = 5.0 * 3.14159265358979323846 / 180.0;
  double min_range_m = 1.0;  ///< drop detections closer than this (DC / leakage)
  /// Minimum SNR (vs. the global noise floor) of each angle peak. Rejects
  /// sidelobe peaks that noise lifts over the secondary-peak threshold in weak
  /// Doppler-leakage cells next to strong static reflectors.
  double min_snr_db = 10.0;
  /// RCS = SNR + 40 log10(R) - rcs_offset_db. Calibrate with a reflector of
  /// known RCS (see calibrate_rcs_offset).
  double rcs_offset_db = 0.0;
};

/// Offset that maps (SNR, range) of a reflector with known RCS onto that RCS.
inline double calibrate_rcs_offset(double measured_snr_db, double range_m, double known_rcs_dbsm) {
  return measured_snr_db + 40.0 * std::log10(range_m) - known_rcs_dbsm;
}

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
  double noise_floor = 0.0;    ///< mean noise power per cell (global estimate)
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
  void set_rcs_offset_db(double offset) { cfg_.rcs_offset_db = offset; }

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
  std::vector<AngleEstimate> angles_;
  struct Candidate {
    std::size_t cell;  ///< index into cells_
    double azimuth_rad;
    double power;
  };
  std::vector<Candidate> candidates_;
  std::vector<std::size_t> cell_first_;  ///< candidates of cell i: [first[i], first[i+1])
  std::vector<char> keep_;
  Grid2D<int> cell_lookup_;
  std::vector<cfloat> steering_;  ///< conjugate steering vector of the current peak  ///< (range,
                                  ///< doppler) -> CFAR cell index, -1 if none
  std::vector<float> scratch_;
};

/// Mean noise power per cell from the median of the power map. The median is
/// robust to targets and clutter; for Gamma(K) noise median ~ K (1 - 1/(9K))^3.
double estimate_noise_floor(const Grid2D<float>& power, std::size_t num_looks,
                            std::vector<float>& scratch);

}  // namespace radar
