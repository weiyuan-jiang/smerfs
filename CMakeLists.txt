cmake_minimum_required(VERSION 3.18)

# ---------------------------------------------------------------------------
# Compilers — must be set before project().
# Load the nag-stack environment module first:
#   ml nag-stack
# That exports FC=nagfor and CC=cc (Apple clang), which CMake picks up here.
# ---------------------------------------------------------------------------
if(NOT DEFINED CMAKE_Fortran_COMPILER AND DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER "$ENV{FC}" CACHE STRING "Fortran compiler")
endif()
if(NOT DEFINED CMAKE_C_COMPILER AND DEFINED ENV{CC})
  set(CMAKE_C_COMPILER "$ENV{CC}" CACHE STRING "C compiler")
endif()

# ---------------------------------------------------------------------------
# Default install prefix: <build_dir>/install
# Set before project() so it wins over CMake's platform default (/usr/local).
# The user can still override with -DCMAKE_INSTALL_PREFIX=/your/path because
# that flag is processed before this file runs.
# ---------------------------------------------------------------------------
if(NOT DEFINED CMAKE_INSTALL_PREFIX OR CMAKE_INSTALL_PREFIX STREQUAL "")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install"
      CACHE PATH "Install prefix (default: <build_dir>/install)" FORCE)
endif()

project(smerfs VERSION 1.0 LANGUAGES C Fortran)

# ---------------------------------------------------------------------------
# Standards
# ---------------------------------------------------------------------------
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Build type
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Compiler flags (applied project-wide)
# ---------------------------------------------------------------------------
# C: optimisation + warnings
add_compile_options(
  $<$<COMPILE_LANGUAGE:C>:-O3>
  $<$<COMPILE_LANGUAGE:C>:-Wall>
)

# Fortran (nagfor):
#   -f2018     accept Fortran 2018 (superset of F90; needed for iso_c_binding)
#   -O3        optimisation
#   -ieee=full do not abort on transient NaN/Inf from C complex arithmetic
add_compile_options(
  $<$<COMPILE_LANGUAGE:Fortran>:-f2018>
  $<$<COMPILE_LANGUAGE:Fortran>:-O3>
  $<$<COMPILE_LANGUAGE:Fortran>:-ieee=full>
)

# ---------------------------------------------------------------------------
# Fortran module output directory (shared by all subdirectories)
# ---------------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# ---------------------------------------------------------------------------
# Installation layout (GNU standard paths)
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# Subdirectories
#   src/        — C kernels (smerfs_c) + Fortran interface (smerfs_f)
#   examples_f/ — Fortran example executable (example1_f)
# ---------------------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(examples_f)

# ---------------------------------------------------------------------------
# Install the Fortran module file
# (targets themselves are installed inside src/CMakeLists.txt)
# ---------------------------------------------------------------------------
install(DIRECTORY "${CMAKE_BINARY_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/smerfs"
  FILES_MATCHING PATTERN "*.mod"
)

# ---------------------------------------------------------------------------
# Install the example executable
# ---------------------------------------------------------------------------
install(TARGETS example1_f
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "============================================================")
message(STATUS "  SMERFS")
message(STATUS "  Version      : ${PROJECT_VERSION}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C compiler   : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}")
message(STATUS "  Fortran      : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}")
message(STATUS "  Install      : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "    libs       -> ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
message(STATUS "    headers    -> ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/smerfs")
message(STATUS "    binaries   -> ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}")
message(STATUS "============================================================")
