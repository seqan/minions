# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/matanatmammadli/Desktop/Bachelorarbeit/minions

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug

# Include any dependencies generated for this target.
include test/api/CMakeFiles/modmer_hash_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/api/CMakeFiles/modmer_hash_test.dir/compiler_depend.make

# Include the progress variables for this target.
include test/api/CMakeFiles/modmer_hash_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/api/CMakeFiles/modmer_hash_test.dir/flags.make

test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o: test/api/CMakeFiles/modmer_hash_test.dir/flags.make
test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o: ../test/api/modmer_hash_test.cpp
test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o: test/api/CMakeFiles/modmer_hash_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o"
	cd /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o -MF CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o.d -o CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o -c /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/api/modmer_hash_test.cpp

test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.i"
	cd /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/api/modmer_hash_test.cpp > CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.i

test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.s"
	cd /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/api/modmer_hash_test.cpp -o CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.s

# Object files for target modmer_hash_test
modmer_hash_test_OBJECTS = \
"CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o"

# External object files for target modmer_hash_test
modmer_hash_test_EXTERNAL_OBJECTS =

test/api/modmer_hash_test: test/api/CMakeFiles/modmer_hash_test.dir/modmer_hash_test.cpp.o
test/api/modmer_hash_test: test/api/CMakeFiles/modmer_hash_test.dir/build.make
test/api/modmer_hash_test: lib/libminions_lib.a
test/api/modmer_hash_test: lib/libgtestd.a
test/api/modmer_hash_test: lib/libgtest_maind.a
test/api/modmer_hash_test: /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/lib/libz.tbd
test/api/modmer_hash_test: /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/lib/libbz2.tbd
test/api/modmer_hash_test: lib/libstrobemer_lib.a
test/api/modmer_hash_test: lib/libgtestd.a
test/api/modmer_hash_test: test/api/CMakeFiles/modmer_hash_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable modmer_hash_test"
	cd /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/modmer_hash_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/api/CMakeFiles/modmer_hash_test.dir/build: test/api/modmer_hash_test
.PHONY : test/api/CMakeFiles/modmer_hash_test.dir/build

test/api/CMakeFiles/modmer_hash_test.dir/clean:
	cd /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api && $(CMAKE_COMMAND) -P CMakeFiles/modmer_hash_test.dir/cmake_clean.cmake
.PHONY : test/api/CMakeFiles/modmer_hash_test.dir/clean

test/api/CMakeFiles/modmer_hash_test.dir/depend:
	cd /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/matanatmammadli/Desktop/Bachelorarbeit/minions /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/api /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/api/CMakeFiles/modmer_hash_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/api/CMakeFiles/modmer_hash_test.dir/depend

