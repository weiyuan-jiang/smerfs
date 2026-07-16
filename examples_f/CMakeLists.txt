# examples_f/CMakeLists.txt
# Included via add_subdirectory() from the root CMakeLists.txt.
# The smerfs_c and smerfs_f targets are defined in src/CMakeLists.txt,
# which is processed first by the root.

# ---------------------------------------------------------------------------
# Executable: example1_f
# Fortran translation of examples/example1.py
# ---------------------------------------------------------------------------
add_executable(example1_f example1_f.F90)

# Output filename: example1_f.x
set_target_properties(example1_f PROPERTIES OUTPUT_NAME example1_f.x)

# smerfs_f carries the module directory and smerfs_c as transitive dependencies
target_link_libraries(example1_f PRIVATE smerfs_f m)

# example1_f.x takes no input files; all parameters are compiled in.
# It writes realisation.csv and correl.csv to the working directory.
