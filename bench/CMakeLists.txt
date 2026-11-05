add_executable(radar_bench bench_pipeline.cpp)
target_link_libraries(radar_bench PRIVATE radar::core benchmark::benchmark)
