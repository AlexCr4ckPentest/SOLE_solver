# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/alex/SOLE_solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alex/SOLE_solver/build

# Include any dependencies generated for this target.
include CMakeFiles/SOLE-solver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SOLE-solver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SOLE-solver.dir/flags.make

CMakeFiles/SOLE-solver.dir/src/main.cpp.o: CMakeFiles/SOLE-solver.dir/flags.make
CMakeFiles/SOLE-solver.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alex/SOLE_solver/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SOLE-solver.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/SOLE-solver.dir/src/main.cpp.o -c /home/alex/SOLE_solver/src/main.cpp

CMakeFiles/SOLE-solver.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SOLE-solver.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alex/SOLE_solver/src/main.cpp > CMakeFiles/SOLE-solver.dir/src/main.cpp.i

CMakeFiles/SOLE-solver.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SOLE-solver.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alex/SOLE_solver/src/main.cpp -o CMakeFiles/SOLE-solver.dir/src/main.cpp.s

# Object files for target SOLE-solver
SOLE__solver_OBJECTS = \
"CMakeFiles/SOLE-solver.dir/src/main.cpp.o"

# External object files for target SOLE-solver
SOLE__solver_EXTERNAL_OBJECTS =

SOLE-solver: CMakeFiles/SOLE-solver.dir/src/main.cpp.o
SOLE-solver: CMakeFiles/SOLE-solver.dir/build.make
SOLE-solver: CMakeFiles/SOLE-solver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/alex/SOLE_solver/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SOLE-solver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SOLE-solver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SOLE-solver.dir/build: SOLE-solver

.PHONY : CMakeFiles/SOLE-solver.dir/build

CMakeFiles/SOLE-solver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SOLE-solver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SOLE-solver.dir/clean

CMakeFiles/SOLE-solver.dir/depend:
	cd /home/alex/SOLE_solver/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alex/SOLE_solver /home/alex/SOLE_solver /home/alex/SOLE_solver/build /home/alex/SOLE_solver/build /home/alex/SOLE_solver/build/CMakeFiles/SOLE-solver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SOLE-solver.dir/depend

