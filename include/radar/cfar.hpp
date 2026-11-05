#pragma once

#include <cstddef>
#include <vector>

#include "radar/types.hpp"

namespace radar {

struct CfarConfig {
  std::size_t train_range = 8;    ///< training cells on each side, range axis
  std::size_t train_doppler = 4;  ///< training cells on each side, Doppler axis
  std::size_t guard_range = 2;
  std::size_t guard_doppler = 2;
  double pfa = 1e-6;  ///< design probability of false alarm per cell
  /// Number of non-coherently summed looks per cell (e.g. receive channels).
  /// Noise in each cell is then Gamma(num_looks)-distributed.
  std::size_t num_looks = 1;
};

struct CfarCell {
  std::size_t range_bin = 0;
  std::size_t doppler_bin = 0;
  float power = 0.0f;
  float noise = 0.0f;  ///< mean training-cell power (local noise estimate)
};

/// 2D cell-averaging CFAR on a power map indexed (range, doppler).
///
/// - Window sums come from a summed-area table: O(1) per cell regardless of
///   window size.
/// - The Doppler axis is circular (velocity aliases), so cells near +/- v_max
///   use wrapped training cells instead of being skipped.
/// - Near range edges the window is clipped and the threshold is recomputed
///   for the reduced number of training cells, so no cell is left untested.
/// - The threshold multiplier is derived from the requested Pfa exactly for
///   square-law detection with `num_looks` integrated looks.
class CaCfar2D {
 public:
  explicit CaCfar2D(const CfarConfig& cfg);

  void detect(const Grid2D<float>& power, std::vector<CfarCell>& out);

  /// Multiplier t such that P(cell > t * sum_of_training_cells) = pfa, for
  /// `num_train` training cells of Gamma(`num_looks`) noise.
  static double threshold_multiplier(std::size_t num_train, std::size_t num_looks, double pfa);

  /// Analytic false-alarm probability for a given multiplier (inverse of the above).
  static double false_alarm_probability(double multiplier, std::size_t num_train,
                                        std::size_t num_looks);

  const CfarConfig& config() const { return cfg_; }

 private:
  double multiplier_for(std::size_t num_train);

  CfarConfig cfg_;
  Grid2D<double> sat_;
  std::vector<double> multiplier_cache_;
};

/// Keeps only cells that are the maximum of their 3x3 neighbourhood
/// (Doppler wraps, range clips), collapsing each target to one cell.
void suppress_non_maxima(const Grid2D<float>& power, const std::vector<CfarCell>& cells,
                         std::vector<CfarCell>& peaks);

}  // namespace radar
