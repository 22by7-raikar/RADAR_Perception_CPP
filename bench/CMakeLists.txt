add_executable(radar_bench bench_pipeline.cpp)
target_link_libraries(radar_bench PRIVATE radar::core rp_options benchmark::benchmark)
