cmake_minimum_required(VERSION 3.18)

# ---------------------------------------------------------------------------
# Compilers must be set before project() is called.
# When built with the nag-stack environment module loaded, the environment
# already exports FC=nagfor and CC=cc (Apple clang).  CMake picks these up
# automatically; the lines below are a safety fallback.
# ---------------------------------------------------------------------------
if(NOT DEFINED CMAKE_Fortran_COMPILER AND DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER "$ENV{FC}" CACHE STRING "Fortran compiler")
endif()
if(NOT DEFINED CMAKE_C_COMPILER AND DEFINED ENV{CC})
  set(CMAKE_C_COMPILER "$ENV{CC}" CACHE STRING "C compiler")
endif()

project(smerfs LANGUAGES C Fortran)

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
# Compiler flags
# ---------------------------------------------------------------------------
# C (Apple clang / gcc)
add_compile_options(
  $<$<COMPILE_LANGUAGE:C>:-O3>
  $<$<COMPILE_LANGUAGE:C>:-Wall>
)

# Fortran (nagfor)
# -f2018   : accept Fortran 2018 (superset of F90; needed for iso_c_binding)
# -C=all   : optional runtime checks (comment out for production)
# -O3      : optimisation
add_compile_options(
  $<$<COMPILE_LANGUAGE:Fortran>:-f2018>
  $<$<COMPILE_LANGUAGE:Fortran>:-O3>
  $<$<COMPILE_LANGUAGE:Fortran>:-ieee=full>
)

# Put all .mod files in one known directory so dependents can find them
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

# ---------------------------------------------------------------------------
# C library  (smerfs_c)
# Sources: cov.c  linalg.c  smerfs.c  ziggurat.c
# zig.h is included by ziggurat.c and lives in the same directory.
# ---------------------------------------------------------------------------
add_library(smerfs_c STATIC
  cov.c
  linalg.c
  smerfs.c
  ziggurat.c
)

target_include_directories(smerfs_c PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

# ---------------------------------------------------------------------------
# Fortran interface library  (smerfs_f)
# Source: smerfs_interface.f90
# ---------------------------------------------------------------------------
add_library(smerfs_f STATIC
  smerfs_interface.f90
)

# Consumers of smerfs_f can USE smerfs_interface by pointing to the module dir
target_include_directories(smerfs_f PUBLIC
  "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

# Link the Fortran interface against the C kernels
target_link_libraries(smerfs_f PUBLIC smerfs_c)

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS smerfs_c smerfs_f
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/smerfs"
  FILES_MATCHING PATTERN "*.mod"
)

install(FILES zig.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/smerfs"
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "------------------------------------------------------------")
message(STATUS "  SMERFS build configuration")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C compiler   : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}")
message(STATUS "  Fortran      : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER}")
message(STATUS "  Module dir   : ${CMAKE_Fortran_MODULE_DIRECTORY}")
message(STATUS "------------------------------------------------------------")
