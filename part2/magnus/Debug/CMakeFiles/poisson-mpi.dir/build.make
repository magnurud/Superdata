# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/magnus/SD/Superdata/part2/magnus

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/magnus/SD/Superdata/part2/magnus/Debug

# Include any dependencies generated for this target.
include CMakeFiles/poisson-mpi.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/poisson-mpi.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/poisson-mpi.dir/flags.make

CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o: CMakeFiles/poisson-mpi.dir/flags.make
CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o: ../poisson-mpi.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/magnus/SD/Superdata/part2/magnus/Debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o   -c /home/magnus/SD/Superdata/part2/magnus/poisson-mpi.c

CMakeFiles/poisson-mpi.dir/poisson-mpi.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/poisson-mpi.dir/poisson-mpi.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/magnus/SD/Superdata/part2/magnus/poisson-mpi.c > CMakeFiles/poisson-mpi.dir/poisson-mpi.c.i

CMakeFiles/poisson-mpi.dir/poisson-mpi.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/poisson-mpi.dir/poisson-mpi.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/magnus/SD/Superdata/part2/magnus/poisson-mpi.c -o CMakeFiles/poisson-mpi.dir/poisson-mpi.c.s

CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.requires:
.PHONY : CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.requires

CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.provides: CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.requires
	$(MAKE) -f CMakeFiles/poisson-mpi.dir/build.make CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.provides.build
.PHONY : CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.provides

CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.provides.build: CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o

# Object files for target poisson-mpi
poisson__mpi_OBJECTS = \
"CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o"

# External object files for target poisson-mpi
poisson__mpi_EXTERNAL_OBJECTS =

poisson-mpi: CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o
poisson-mpi: CMakeFiles/poisson-mpi.dir/build.make
poisson-mpi: libcommon.a
poisson-mpi: /usr/lib/libmpi.so
poisson-mpi: /usr/lib/x86_64-linux-gnu/libdl.so
poisson-mpi: /usr/lib/x86_64-linux-gnu/libhwloc.so
poisson-mpi: CMakeFiles/poisson-mpi.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable poisson-mpi"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poisson-mpi.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/poisson-mpi.dir/build: poisson-mpi
.PHONY : CMakeFiles/poisson-mpi.dir/build

CMakeFiles/poisson-mpi.dir/requires: CMakeFiles/poisson-mpi.dir/poisson-mpi.c.o.requires
.PHONY : CMakeFiles/poisson-mpi.dir/requires

CMakeFiles/poisson-mpi.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/poisson-mpi.dir/cmake_clean.cmake
.PHONY : CMakeFiles/poisson-mpi.dir/clean

CMakeFiles/poisson-mpi.dir/depend:
	cd /home/magnus/SD/Superdata/part2/magnus/Debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/magnus/SD/Superdata/part2/magnus /home/magnus/SD/Superdata/part2/magnus /home/magnus/SD/Superdata/part2/magnus/Debug /home/magnus/SD/Superdata/part2/magnus/Debug /home/magnus/SD/Superdata/part2/magnus/Debug/CMakeFiles/poisson-mpi.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/poisson-mpi.dir/depend

