#include <benchmark/benchmark.h>

#include "radar/pipeline.hpp"
#include "radar/simulator.hpp"

namespace {

using namespace radar;

std::vector<Target> scene() {
  return {{12.0, 0.0, -0.4, 15.0},   {35.0, -4.5, 0.17, 0.0},   {60.0, 12.0, -0.14, 20.0},
          {60.0, -12.0, 0.26, 20.0}, {95.0, -25.0, 0.52, 10.0}, {150.0, 8.0, -0.7, 25.0}};
}

const DataCube& frame() {
  static const DataCube f = [] {
    Simulator sim{RadarConfig{}};
    return sim.simulate_frame(scene());
  }();
  return f;
}

void BM_RangeDopplerFft(benchmark::State& state) {
  RadarConfig cfg;
  RangeDopplerProcessor rd(cfg);
  DataCube spec(cfg.num_rx, cfg.chirps_per_frame, cfg.samples_per_chirp);
  const DataCube& adc = frame();  // simulate outside the timed loop
  for (auto _ : state) {
    rd.process(adc, spec);
    benchmark::DoNotOptimize(spec.data());
  }
  state.SetItemsProcessed(state.iterations() * static_cast<long>(frame().size()));
}
BENCHMARK(BM_RangeDopplerFft)->Unit(benchmark::kMicrosecond);

void BM_NonCoherentIntegration(benchmark::State& state) {
  RadarConfig cfg;
  RangeDopplerProcessor rd(cfg);
  DataCube spec;
  rd.process(frame(), spec);
  Grid2D<float> power;
  for (auto _ : state) {
    integrate_power(spec, power);
    benchmark::DoNotOptimize(power.data());
  }
}
BENCHMARK(BM_NonCoherentIntegration)->Unit(benchmark::kMicrosecond);

void BM_CaCfar2D(benchmark::State& state) {
  RadarConfig cfg;
  RangeDopplerProcessor rd(cfg);
  DataCube spec;
  rd.process(frame(), spec);
  Grid2D<float> power;
  integrate_power(spec, power);
  CfarConfig cc;
  cc.num_looks = cfg.num_rx;
  CaCfar2D cfar(cc);
  std::vector<CfarCell> cells;
  for (auto _ : state) {
    cfar.detect(power, cells);
    benchmark::DoNotOptimize(cells.data());
  }
  state.SetItemsProcessed(state.iterations() * static_cast<long>(power.size()));
}
BENCHMARK(BM_CaCfar2D)->Unit(benchmark::kMicrosecond);

void BM_FullPipeline(benchmark::State& state) {
  RadarPipeline pipe{RadarConfig{}};
  FrameResult res;
  const DataCube& adc = frame();
  for (auto _ : state) {
    pipe.process(adc, res);
    benchmark::DoNotOptimize(res.detections.data());
  }
}
BENCHMARK(BM_FullPipeline)->Unit(benchmark::kMicrosecond);

}  // namespace

BENCHMARK_MAIN();
