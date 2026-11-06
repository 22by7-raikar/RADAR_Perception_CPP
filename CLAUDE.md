# Working on this repository

FMCW radar perception in C++17: raw IQ frames → range-Doppler → CFAR → angle
→ point cloud → tracking, plus ego-motion, radar-camera fusion, TDM-MIMO 4D,
real-world impairments, ROS 2 wrappers. Everything runs on a physically based
simulator; there is no real dataset yet. README.md is the full technical
write-up; read the section for whatever you touch.

## Build and verify

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DRP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build -L unit --output-on-failure          # ~100 tests, seconds
(cd build/tests && ./radar_regression)                       # 42 metrics vs thresholds, ~1 min
./build/bench/radar_bench                                    # Google Benchmark
clang-format-18 --dry-run --Werror $(git ls-files '*.cpp' '*.hpp')
```

CI (`.github/workflows/ci.yml`) runs: clang-format, GCC and Clang Release
with `-Werror` (unit + regression + demos + benchmark smoke), GCC
ASan/UBSan, aarch64 cross-compile with tests under QEMU, and ROS 2 Jazzy
(4D array) / Humble (8 Rx) launch smoke tests. A change is done when all of
that would pass; run the Clang build too if you can, it has caught real
issues (unused lambda captures, different noise realisations in borderline
tests).

## Layout

- `include/radar/`, `src/`: the library (`radar::core`). Signal chain:
  `config` (incl. TDM-MIMO layout), `simulator`, `range_doppler`
  (pocketfft), `cfar` (CA), `angle` (virtual array, elevation), `pipeline`
  (owns all buffers; orchestrates everything per frame). Objects:
  `clustering`, `tracker`, `assignment`. Ego: `ego_velocity`, `imu`,
  `ego_estimator`. Camera: `camera`, `calibration`, `fusion`. Impairments:
  `impairments` (interference, clutter, multipath models), `mitigation`
  (blanking + IMAT, OS-CFAR, guardrail ghost filter),
  `impairment_experiments` (shared by demo and regression).
- `apps/*_demo.cpp` write CSVs; `tools/plot_*.py` turn them into
  `docs/images/*.png` (keep the existing palette/rcParams).
- `tests/` unit tests (GoogleTest, label `unit`); `tests/regression/`
  fixed-seed scenarios checked against `thresholds.txt` (bound + measured
  value in a comment; GCC and Clang can differ slightly, bounds cover both).
- `ros2/`: thin nodes over the library; built against an installed core
  (`cmake --install`, `find_package(radar_perception)`).

## Conventions

- C++17, clang-format (Google-based, 100 cols), warnings clean under GCC and
  Clang with `-Wall -Wextra -Wpedantic -Wshadow`.
- No allocation in the steady-state per-frame path: buffers live in the
  owning class and are reused. Timing goes in `StageTimings`.
- Statistical claims need a test: thresholds derived for a Pfa get a Monte
  Carlo false-alarm test; estimators with an uncertainty get a
  predicted-vs-measured spread test.
- New behaviour gets unit tests and, if it changes accuracy, a regression
  metric with a bound. Never loosen a bound without saying why in the
  comment next to it.
- Docs: README sections state measured numbers, what failed along the way
  and what is still not solved. Keep them about the engineering only; no
  references to specific companies or job postings.
- Git: commit with the repository's configured identity only; no
  co-author trailers. Don't push without being asked.

## Next task: CUDA backend (RTX 3060 laptop, compute capability 8.6)

Goal: GPU implementations of the heavy, regular stages, verified against
the CPU path and benchmarked honestly (including transfers).

Plan, in order:

1. **Build plumbing.** `RP_ENABLE_CUDA` option (default OFF so CI and the
   ARM build are unaffected), `enable_language(CUDA)`,
   `find_package(CUDAToolkit)`, `CMAKE_CUDA_ARCHITECTURES=86`, a separate
   `radar_cuda` target linking `radar::core` and `CUDA::cufft`. GPU tests get
   their own label (`gpu`) and only build when the option is on.
2. **Range-Doppler on GPU.** Window kernel (range × Doppler taper incl. the
   (-1)^m Doppler centring, see `RangeDopplerProcessor`), batched cuFFT 2D
   C2C over all channels, then a non-coherent integration kernel → power
   map. Keep the cube layout `[channel][chirp][sample]`.
3. **CA-CFAR on GPU.** Either a summed-area table (scan) or a direct
   window-sum kernel with shared-memory tiles; Doppler wraps, range clips,
   threshold multipliers precomputed on the host per training-set size
   (`CaCfar2D::threshold_multiplier`). Output: compacted cell list.
4. **Correctness tests.** GPU vs CPU on simulated frames: spectra within
   float tolerance, power maps within ~1e-4 relative, CFAR cell sets equal
   except for cells within a hair of the threshold (count those, assert
   tiny).
5. **Benchmarks.** CPU vs GPU per stage with cudaEvents, with and without
   host↔device copies, pinned host memory, a stream to overlap the copy of
   frame N+1 with compute of frame N. Report three sizes: the default
   (8 × 128 × 256, expect the GPU *not* to win once copies count), the
   3 × 4 TDM-MIMO, and a large imaging-radar config (e.g. 12 Tx × 16 Rx
   = 192 virtual channels, 256 chirps, 512 samples; needs a new
   `RadarConfig` factory whose virtual azimuth row is uniform, see
   `analyse_virtual_array`). That large case is where the GPU should pay off.
6. **Optional:** interference mitigation (per-chirp FFTs + IMAT) is the
   most expensive CPU stage when interference is present; it batches well
   on the GPU.
7. **README:** a "GPU" section with the measured table, the GPU model,
   driver/CUDA version, and the honest crossover point.

Environment notes: Linux or WSL2 (Ubuntu 22.04/24.04) with the CUDA
toolkit 12.x is the smoothest; native Windows needs MSVC and the CMake
presets adjusted. Check `nvidia-smi` and `nvcc --version` first.
