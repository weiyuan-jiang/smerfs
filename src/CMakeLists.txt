esma_set_this(OVERRIDE smerfs_f)

set (srcs
  cov.c
  linalg.c
  smerfs.c
  ziggurat.c
  smerfs_fenv.c      # FP-environment wrappers (feholdexcept/fesetenv)
  smerfs_interface.f90
)

esma_add_library (${this}
  SOURCES ${srcs})

# Expose the module directory to consumers of this target
target_include_directories(smerfs_f PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/include"
)
