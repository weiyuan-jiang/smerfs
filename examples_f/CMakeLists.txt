cmake_minimum_required(VERSION 3.18)

# ---------------------------------------------------------------------------
# Compilers: pick up FC/CC set by "ml nag-stack" before calling cmake
# ---------------------------------------------------------------------------
if(NOT DEFINED CMAKE_Fortran_COMPILER AND DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER "$ENV{FC}" CACHE STRING "Fortran compiler")
endif()
if(NOT DEFINED CMAKE_C_COMPILER AND DEFINED ENV{CC})
  set(CMAKE_C_COMPILER "$ENV{CC}" CACHE STRING "C compiler")
endif()

project(smerfs_examples LANGUAGES C Fortran)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Compiler flags (set per-language on the target to avoid shell quoting issues)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Point at the src/ build tree to get libsmerfs_c, libsmerfs_f and the .mod
# ---------------------------------------------------------------------------
# Allow the user to override:  cmake -DSMERFS_SRC_BUILD=<path>
if(NOT DEFINED SMERFS_SRC_BUILD)
  # Default: assume "cmake -S smerfs/examples_f -B build/examples_f"
  # and the src libraries were built in build/src (sibling directory)
  get_filename_component(_default_src_build
    "${CMAKE_CURRENT_SOURCE_DIR}/../build/src" ABSOLUTE)
  set(SMERFS_SRC_BUILD "${_default_src_build}"
      CACHE PATH "Path to the smerfs/src build directory")
endif()

message(STATUS "Looking for smerfs libraries in: ${SMERFS_SRC_BUILD}")

# Import the pre-built static libraries
find_library(SMERFS_F_LIB   NAMES smerfs_f   HINTS "${SMERFS_SRC_BUILD}" REQUIRED)
find_library(SMERFS_C_LIB   NAMES smerfs_c   HINTS "${SMERFS_SRC_BUILD}" REQUIRED)

# Module files
set(SMERFS_MOD_DIR "${SMERFS_SRC_BUILD}/modules")

message(STATUS "  smerfs_f : ${SMERFS_F_LIB}")
message(STATUS "  smerfs_c : ${SMERFS_C_LIB}")
message(STATUS "  modules  : ${SMERFS_MOD_DIR}")

# ---------------------------------------------------------------------------
# Executable: example1_f
# ---------------------------------------------------------------------------
add_executable(example1_f example1_f.F90)

target_include_directories(example1_f PRIVATE "${SMERFS_MOD_DIR}")

target_compile_options(example1_f PRIVATE
  # -ieee=full: do not abort on FPE from C library calls (transient NaN/Inf
  #             in complex hypergeometric series are benign and resolve)
  $<$<COMPILE_LANGUAGE:Fortran>:-f2018>
  $<$<COMPILE_LANGUAGE:Fortran>:-O3>
  $<$<COMPILE_LANGUAGE:Fortran>:-ieee=full>
  $<$<COMPILE_LANGUAGE:C>:-O3>
  $<$<COMPILE_LANGUAGE:C>:-Wall>
)

target_link_libraries(example1_f PRIVATE
  "${SMERFS_F_LIB}"
  "${SMERFS_C_LIB}"
  m   # libm for sin/cos/sqrt in C kernels
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "------------------------------------------------------------")
message(STATUS "  smerfs examples_f build")
message(STATUS "  C compiler   : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}")
message(STATUS "  Fortran      : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}")
message(STATUS "------------------------------------------------------------")
