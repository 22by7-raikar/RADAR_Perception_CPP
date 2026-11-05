# Cross-compile for 64-bit ARM Linux (Jetson / Cortex-A class SoCs).
#
#   cmake -S . -B build-arm64 -G Ninja \
#     -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/aarch64-linux-gnu.cmake \
#     -DRP_TARGET_CPU=cortex-a78ae          # optional, e.g. Jetson Orin
#
# Tests run on the build host through QEMU user-mode emulation when available.
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR aarch64)

set(CMAKE_C_COMPILER aarch64-linux-gnu-gcc)
set(CMAKE_CXX_COMPILER aarch64-linux-gnu-g++)

set(CMAKE_FIND_ROOT_PATH /usr/aarch64-linux-gnu)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
# Header-only dependencies (Eigen) may come from the host.
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH)

find_program(QEMU_AARCH64 qemu-aarch64)
if(QEMU_AARCH64)
  set(CMAKE_CROSSCOMPILING_EMULATOR ${QEMU_AARCH64} -L /usr/aarch64-linux-gnu)
endif()
