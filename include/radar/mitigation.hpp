#pragma once

// Counter-measures for the impairments in impairments.hpp.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "radar/cfar.hpp"
#include "radar/config.hpp"
#include "radar/geometry.hpp"
#include "radar/types.hpp"

namespace radar {

// --- Interference ------------------------------------------------------------

struct InterferenceMitigationConfig {
  bool enabled = false;
  /// Spectral bins this far above the chirp's median spectrum are treated as
  /// target echoes (narrowband) and notched out before burst detection.
  double notch_db = 10.0;
  /// A sample of the notched residual is interference when its power exceeds
  /// the residual's median by this much (single channel: noise alone exceeds
  /// 13 dB over the median with probability ~1e-6).
  double threshold_db = 13.0;
  std::size_t guard_samples = 2;  ///< also blank this many samples either side
  /// Refill blanked samples instead of leaving zeros (0 = zeroing only).
  /// Zeroing cuts gaps into strong echoes at a different place in every
  /// chirp; that modulation smears them over all Doppler bins and nearby
  /// ranges, where CFAR detects the smear. Reconstruction is iterative
  /// thresholding (IMAT), with thresholds evenly spaced in dB from the
  /// chirp's spectral peak down to the notch level.
  std::size_t reconstruct_iterations = 5;
};

/// Detects interference bursts in the time-domain ADC samples and zeroes
/// them, chirp by chirp and before the range FFT. Returns the number of
/// (chirp, sample) slots blanked.
///
/// Thresholding the raw samples fails when a strong near target dominates
/// them (its echo can exceed the interference). Echoes are tones and bursts
/// are broadband, so detection runs on a residual: the chirp's spectrum with
/// strong bins notched out, transformed back to time. Blanking itself is
/// applied to the original samples on every receiver of that chirp, and the
/// gap is refilled from the strong tones (see `reconstruct_iterations`).
std::size_t mitigate_interference(const RadarConfig& radar, const InterferenceMitigationConfig& cfg,
                                  DataCube& adc, std::vector<float>& scratch);

/// The detector's input for one chirp: Hann-windowed samples with spectral
/// bins `notch_db` above the median removed, back in time, as power.
std::vector<float> notched_residual(const cfloat* samples, std::size_t n, double notch_db);

// --- Ordered-statistic CFAR --------------------------------------------------

/// 2D OS-CFAR: the noise estimate is the k-th smallest training cell
/// (k = rank * N) instead of the mean. A strong target or a clutter edge in
/// the training window only shifts the order statistic by a few places, so
/// a weak target next to it is not masked, at the cost of a sort per cell
/// and a small CFAR loss in homogeneous noise.
///
/// The threshold multiplier is solved for the requested Pfa exactly
/// (numerical integration over the order-statistic density of Gamma(K)
/// noise). Only cells above kPrescreen x the homogeneous-noise threshold,
/// judged against the map median, are tested; the cells skipped this way
/// would need an order statistic far below the median to pass. That holds
/// while the noise floor is roughly flat across the map.
class OsCfar2D {
 public:
  explicit OsCfar2D(const CfarConfig& cfg, double rank = 0.75);

  void detect(const Grid2D<float>& power, std::vector<CfarCell>& out);

  /// Solves the thresholds for every training-set size a map with
  /// `range_bins` rows produces (full and range-clipped windows), so the first
  /// detect() call does not pay for them (tens of ms each).
  void prepare(std::size_t range_bins);

  /// Multiplier t with P(cell > t * k-th smallest of N training cells) = pfa.
  static double threshold_multiplier(std::size_t num_train, std::size_t k, std::size_t num_looks,
                                     double pfa);
  static double false_alarm_probability(double multiplier, std::size_t num_train, std::size_t k,
                                        std::size_t num_looks);

  std::size_t rank_for(std::size_t num_train) const;

  static constexpr double kPrescreen = 0.8;

 private:
  double multiplier_for(std::size_t num_train);

  CfarConfig cfg_;
  double rank_;
  std::vector<double> multiplier_cache_;
  std::vector<float> train_;
  std::vector<float> scratch_;
};

// --- Multipath ghosts --------------------------------------------------------

/// Line y = slope * x + offset in the sensor frame.
struct ReflectorLine {
  bool valid = false;
  double slope = 0.0;
  double offset = 0.0;
  std::size_t inliers = 0;
  /// Signed distance of a point from the line, positive on the sensor's side.
  double sensor_side_distance(double x, double y) const;
};

struct GhostFilterConfig {
  double static_velocity_tol_mps = 0.5;
  double min_lateral_m = 1.5;  ///< reflector candidates must be off the ego path
  double max_slope = 0.25;     ///< roughly parallel to the direction of travel
  double inlier_distance_m = 0.4;
  std::size_t min_inliers = 8;
  std::size_t ransac_iterations = 100;
  double behind_margin_m = 1.0;  ///< how far past the line a ghost must lie
  bool require_mirror_match = true;
  double mirror_gate_m = 3.0;
  double mirror_velocity_gate_mps = 2.0;
  std::uint32_t seed = 1;
};

/// Fits a reflector (guardrail, wall) on one side (`side` = +1 left, -1
/// right) to stationary detections with RANSAC. `ego_velocity` is the
/// sensor's velocity in its own frame, used to tell static returns apart.
ReflectorLine fit_reflector(const std::vector<Detection>& dets, const Vec2& ego_velocity, int side,
                            const GhostFilterConfig& cfg = {});

/// Flags moving detections that lie behind a fitted reflector and, if
/// required, mirror a detection on the sensor's side (a real object whose
/// reflection they are). `ghost` is resized to dets.size().
void flag_multipath_ghosts(const std::vector<Detection>& dets,
                           const std::vector<ReflectorLine>& reflectors, const Vec2& ego_velocity,
                           const GhostFilterConfig& cfg, std::vector<char>& ghost);

}  // namespace radar
