add_executable(radar_tests
  test_angle.cpp
  test_assignment.cpp
  test_cfar.cpp
  test_clustering.cpp
  test_config.cpp
  test_fusion.cpp
  test_ego_motion.cpp
  test_io.cpp
  test_mimo.cpp
  test_pipeline.cpp
  test_range_doppler.cpp
  test_scene.cpp
  test_simulator.cpp
  test_tracker.cpp
  test_window.cpp
)
target_link_libraries(radar_tests PRIVATE radar::core rp_options GTest::gtest_main)
target_compile_definitions(radar_tests PRIVATE RP_SCENARIO_DIR="${PROJECT_SOURCE_DIR}/scenarios")

include(GoogleTest)
gtest_discover_tests(radar_tests DISCOVERY_TIMEOUT 120 PROPERTIES LABELS unit)

# Scenario-level regression suite: fixed seeds, metrics checked against
# tests/regression/thresholds.txt. Slow under sanitizers; select with -L regression.
add_executable(radar_regression regression/test_regression.cpp)
target_link_libraries(radar_regression PRIVATE radar::core rp_options GTest::gtest_main)
target_compile_definitions(radar_regression PRIVATE
  RP_SCENARIO_DIR="${PROJECT_SOURCE_DIR}/scenarios"
  RP_THRESHOLDS_FILE="${CMAKE_CURRENT_SOURCE_DIR}/regression/thresholds.txt")
gtest_discover_tests(radar_regression DISCOVERY_TIMEOUT 120 PROPERTIES LABELS regression TIMEOUT 600)
