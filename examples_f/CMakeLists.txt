# examples_f/CMakeLists.txt
# Included via add_subdirectory() from the root CMakeLists.txt.
# The smerfs_c and smerfs_f targets are defined in src/CMakeLists.txt,
# which is processed first by the root.

# ---------------------------------------------------------------------------
# Executable: example1_f
# Fortran translation of examples/example1.py
# ---------------------------------------------------------------------------
add_executable(example1_f example1_f.F90)

# smerfs_f carries the module directory and smerfs_c as transitive dependencies
target_link_libraries(example1_f PRIVATE smerfs_f m)

# example1_f writes output CSVs to CWD; no extra include paths needed beyond
# what smerfs_f already propagates (module dir via its PUBLIC include).
