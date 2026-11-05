add_executable(radar_tests
  test_angle.cpp
  test_cfar.cpp
  test_config.cpp
  test_io.cpp
  test_pipeline.cpp
  test_range_doppler.cpp
  test_window.cpp
)
target_link_libraries(radar_tests PRIVATE radar::core GTest::gtest_main)
target_compile_definitions(radar_tests PRIVATE RP_SCENARIO_DIR="${PROJECT_SOURCE_DIR}/scenarios")

include(GoogleTest)
gtest_discover_tests(radar_tests)
