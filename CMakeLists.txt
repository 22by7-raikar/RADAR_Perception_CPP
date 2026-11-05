cmake_minimum_required(VERSION 3.20)
project(radar_perception VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(RP_BUILD_TESTS "Build unit tests" ON)
option(RP_BUILD_BENCHMARKS "Build Google Benchmark suite" ON)
option(RP_BUILD_APPS "Build demo applications" ON)
option(RP_ENABLE_SANITIZERS "Build with AddressSanitizer + UBSan" OFF)
option(RP_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(RP_NATIVE_ARCH "Optimise for the build machine (-march=native)" OFF)

include(cmake/Dependencies.cmake)

# Common compile settings, applied to every first-party target.
add_library(rp_options INTERFACE)
target_compile_options(rp_options INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor)
if(RP_WARNINGS_AS_ERRORS)
  target_compile_options(rp_options INTERFACE -Werror)
endif()
if(RP_NATIVE_ARCH)
  target_compile_options(rp_options INTERFACE -march=native)
endif()
if(RP_ENABLE_SANITIZERS)
  target_compile_options(rp_options INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(rp_options INTERFACE -fsanitize=address,undefined)
endif()

add_library(radar_core
  src/angle.cpp
  src/assignment.cpp
  src/cfar.cpp
  src/clustering.cpp
  src/config.cpp
  src/ego_estimator.cpp
  src/ego_velocity.cpp
  src/evaluation.cpp
  src/imu.cpp
  src/io.cpp
  src/pipeline.cpp
  src/range_doppler.cpp
  src/scenarios.cpp
  src/scene.cpp
  src/simulator.cpp
  src/tracker.cpp
  src/tracking_metrics.cpp
  src/window.cpp
)
add_library(radar::core ALIAS radar_core)
target_include_directories(radar_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(radar_core PUBLIC rp_options Eigen3::Eigen PRIVATE pocketfft)

if(RP_BUILD_APPS)
  add_executable(radar_demo apps/radar_demo.cpp)
  target_link_libraries(radar_demo PRIVATE radar::core)
  add_executable(tracking_demo apps/tracking_demo.cpp)
  target_link_libraries(tracking_demo PRIVATE radar::core)
  add_executable(odometry_demo apps/odometry_demo.cpp)
  target_link_libraries(odometry_demo PRIVATE radar::core)
endif()

if(RP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(RP_BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
