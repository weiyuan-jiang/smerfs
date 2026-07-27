# ---------------------------------------------------------------------------
# Executable: example1_f.x
# Fortran translation of examples/example1.py
# ---------------------------------------------------------------------------
set (exe example1_f.x)

ecbuild_add_executable(
   TARGET ${exe}
   SOURCES example1_f.F90)

# smerfs_f carries the module directory and smerfs_c as transitive dependencies
target_link_libraries(${exe} PRIVATE GEOSlandpert_GridComp MAPL m)
target_include_directories (${exe} PUBLIC $<BUILD_INTERFACE:${MAPL_SOURCE_DIR}/include>)
