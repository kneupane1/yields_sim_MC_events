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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build

# Include any dependencies generated for this target.
include CMakeFiles/clas12_mc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/clas12_mc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/clas12_mc.dir/flags.make

CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.o: CMakeFiles/clas12_mc.dir/flags.make
CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.o: ../src/exe/clas12_mc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.o -c /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/src/exe/clas12_mc.cpp

CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/src/exe/clas12_mc.cpp > CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.i

CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/src/exe/clas12_mc.cpp -o CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.s

# Object files for target clas12_mc
clas12_mc_OBJECTS = \
"CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.o"

# External object files for target clas12_mc
clas12_mc_EXTERNAL_OBJECTS =

clas12_mc: CMakeFiles/clas12_mc.dir/src/exe/clas12_mc.cpp.o
clas12_mc: CMakeFiles/clas12_mc.dir/build.make
clas12_mc: libclas12lib.a
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libCore.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libImt.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libRIO.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libNet.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libHist.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libGraf.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libGraf3d.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libGpad.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libROOTDataFrame.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libTree.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libTreePlayer.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libRint.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libPostscript.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libMatrix.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libPhysics.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libMathCore.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libThread.so
clas12_mc: /usr/local/Cellar/root/6.22.08/lib/root/libMultiProc.so
clas12_mc: CMakeFiles/clas12_mc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable clas12_mc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/clas12_mc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/clas12_mc.dir/build: clas12_mc

.PHONY : CMakeFiles/clas12_mc.dir/build

CMakeFiles/clas12_mc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/clas12_mc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/clas12_mc.dir/clean

CMakeFiles/clas12_mc.dir/depend:
	cd /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build /Users/krishnaneupane/Documents/GitHub/yields_sim_MC_events/build/CMakeFiles/clas12_mc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/clas12_mc.dir/depend

