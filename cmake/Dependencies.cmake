include(FetchContent)

# pocketfft: header-only FFT (BSD-3-Clause). Pinned to a commit on the C++ branch.
FetchContent_Declare(pocketfft
  GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
  GIT_TAG        4d865cba24359b05b37da4e75aee24498c3b78e7
  GIT_SHALLOW    FALSE
)
FetchContent_MakeAvailable(pocketfft)
add_library(pocketfft INTERFACE)
target_include_directories(pocketfft SYSTEM INTERFACE ${pocketfft_SOURCE_DIR})
# Threading is managed by the caller; keep FFTs single-threaded and deterministic.
target_compile_definitions(pocketfft INTERFACE POCKETFFT_NO_MULTITHREADING)

if(RP_BUILD_TESTS)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
    GIT_SHALLOW    TRUE
  )
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

if(RP_BUILD_BENCHMARKS)
  FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        v1.9.1
    GIT_SHALLOW    TRUE
  )
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(benchmark)
endif()
