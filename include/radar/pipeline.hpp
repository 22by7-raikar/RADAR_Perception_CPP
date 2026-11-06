#pragma once

#include <cmath>
#include <vector>

#include "radar/angle.hpp"
#include "radar/cfar.hpp"
#include "radar/config.hpp"
#include "radar/mitigation.hpp"
#include "radar/range_doppler.hpp"
#include "radar/types.hpp"

namespace radar {

enum class CfarType { kCellAveraging, kOrderedStatistic };

struct PipelineConfig {
  /// Time-domain blanking of interference bursts before the range FFT.
  InterferenceMitigationConfig interference;
  RangeDopplerConfig range_doppler;
  CfarConfig cfar;  ///< num_looks is overridden with the number of channels
  CfarType cfar_type = CfarType::kCellAveraging;
  double os_cfar_rank = 0.75;  ///< OS-CFAR only: order statistic as a fraction of N
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
  /// TDM-MIMO only. Transmitters fire one after another, so a moving target's
  /// phase advances between them and corrupts the virtual array; undo it
  /// with the cell's Doppler before angle estimation.
  bool tdm_motion_compensation = true;
  /// TDM-MIMO only. The compensation phase differs between velocity
  /// hypotheses v + k * 2 v_max that the Doppler FFT cannot tell apart; pick
  /// the one that focuses the azimuth beam best. Recovers the
  /// single-transmitter velocity span (num_tx times wider).
  bool velocity_unfolding = true;
  /// Arrays with elevation only: drop peaks outside the antenna's elevation
  /// field of view. Real patch arrays barely radiate there, so such estimates
  /// come from corrupted cells (several unresolved reflectors that defeat the
  /// unfolding test), not from real reflectors.
  double max_abs_elevation_rad = 30.0 * 3.14159265358979323846 / 180.0;
  /// Two velocity-dependent range biases, removed with the (unfolded)
  /// velocity: the Doppler shift adds to the beat frequency, so a target
  /// appears v * fc / S further away (0.25 m at 40 m/s with the default
  /// chirp); and the FFT measures the range averaged over the frame, which is
  /// referred back to the frame start (the cube's timestamp).
  bool range_doppler_correction = true;
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
  std::size_t interference_samples_blanked = 0;  ///< (chirp, sample) slots zeroed
  std::size_t cfar_cells = 0;                    ///< cells above threshold before peak grouping
  double noise_floor = 0.0;                      ///< mean noise power per cell (global estimate)
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
  OsCfar2D os_cfar_;
  AngleEstimator angle_;

  DataCube adc_work_;  ///< copy of the frame for interference blanking
  DataCube spectrum_;
  Grid2D<float> power_;
  std::vector<CfarCell> cells_;
  std::vector<CfarCell> peaks_;
  std::vector<AngleEstimate> angles_;
  struct Candidate {
    std::size_t cell;  ///< index into cells_
    double azimuth_rad;
    double elevation_rad;
    double power;
    int fold;  ///< velocity hypothesis k of the cell
  };
  std::vector<Candidate> candidates_;
  std::vector<std::size_t> cell_first_;  ///< candidates of cell i: [first[i], first[i+1])
  std::vector<char> keep_;
  Grid2D<int> cell_lookup_;
  std::vector<cfloat> steering_;  ///< conjugate steering vector of the current peak
  std::vector<cfloat> snapshot_;  ///< channels of the current cell
  std::vector<cfloat> compensated_;
  std::vector<AngleEstimate> trial_angles_;
  std::vector<float> scratch_;

  /// Angle peaks of one cell into angles_; returns the chosen velocity fold.
  int estimate_cell_angles(std::size_t range_bin, std::size_t doppler_bin);
  /// How well hypothesis-compensated channels (in compensated_) fit one
  /// direction; higher is better.
  double fold_score(double peak_power) const;
  /// Per-channel phase that undoes TDM motion for Doppler bin `d`, fold `k`.
  cfloat tdm_compensation(std::size_t channel, std::size_t doppler_bin, int fold) const;
};

/// Mean noise power per cell from the median of the power map. The median is
/// robust to targets and clutter; for Gamma(K) noise median ~ K (1 - 1/(9K))^3.
double estimate_noise_floor(const Grid2D<float>& power, std::size_t num_looks,
                            std::vector<float>& scratch);

}  // namespace radar
