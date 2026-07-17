# src/CMakeLists.txt
# Included via add_subdirectory() from the root CMakeLists.txt.
# Compiler flags, build type, and install paths are all set by the root.

# ---------------------------------------------------------------------------
# C library: smerfs_c
# Sources: cov.c  linalg.c  smerfs.c  ziggurat.c
# zig.h is included by ziggurat.c and lives in this directory.
# ---------------------------------------------------------------------------
add_library(smerfs_c STATIC
  cov.c
  linalg.c
  smerfs.c
  ziggurat.c
  smerfs_fenv.c      # FP-environment wrappers (feholdexcept/fesetenv)
)

target_include_directories(smerfs_c PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smerfs>"
)

# ---------------------------------------------------------------------------
# Fortran interface library: smerfs_f
# Source: smerfs_interface.f90
# ---------------------------------------------------------------------------
add_library(smerfs_f STATIC
  smerfs_interface.f90
)

# Expose the module directory to consumers of this target
target_include_directories(smerfs_f PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smerfs>"
)

# Link the Fortran wrapper against the C kernels
target_link_libraries(smerfs_f PUBLIC smerfs_c)

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------
install(TARGETS smerfs_c smerfs_f
  EXPORT smerfs_targets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

install(FILES zig.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/smerfs"
)
