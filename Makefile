# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jack/Desktop/gitmt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jack/Desktop/gitmt

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jack/Desktop/gitmt/CMakeFiles /home/jack/Desktop/gitmt/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jack/Desktop/gitmt/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named mt

# Build rule for target.
mt: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mt
.PHONY : mt

# fast build rule for target.
mt/fast:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/build
.PHONY : mt/fast

configreader.o: configreader.cpp.o
.PHONY : configreader.o

# target to build an object file
configreader.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/configreader.cpp.o
.PHONY : configreader.cpp.o

configreader.i: configreader.cpp.i
.PHONY : configreader.i

# target to preprocess a source file
configreader.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/configreader.cpp.i
.PHONY : configreader.cpp.i

configreader.s: configreader.cpp.s
.PHONY : configreader.s

# target to generate assembly for a file
configreader.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/configreader.cpp.s
.PHONY : configreader.cpp.s

dcdio.o: dcdio.cpp.o
.PHONY : dcdio.o

# target to build an object file
dcdio.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/dcdio.cpp.o
.PHONY : dcdio.cpp.o

dcdio.i: dcdio.cpp.i
.PHONY : dcdio.i

# target to preprocess a source file
dcdio.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/dcdio.cpp.i
.PHONY : dcdio.cpp.i

dcdio.s: dcdio.cpp.s
.PHONY : dcdio.s

# target to generate assembly for a file
dcdio.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/dcdio.cpp.s
.PHONY : dcdio.cpp.s

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/main.cpp.s
.PHONY : main.cpp.s

num_test.o: num_test.cpp.o
.PHONY : num_test.o

# target to build an object file
num_test.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/num_test.cpp.o
.PHONY : num_test.cpp.o

num_test.i: num_test.cpp.i
.PHONY : num_test.i

# target to preprocess a source file
num_test.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/num_test.cpp.i
.PHONY : num_test.cpp.i

num_test.s: num_test.cpp.s
.PHONY : num_test.s

# target to generate assembly for a file
num_test.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/num_test.cpp.s
.PHONY : num_test.cpp.s

parameters.o: parameters.cpp.o
.PHONY : parameters.o

# target to build an object file
parameters.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/parameters.cpp.o
.PHONY : parameters.cpp.o

parameters.i: parameters.cpp.i
.PHONY : parameters.i

# target to preprocess a source file
parameters.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/parameters.cpp.i
.PHONY : parameters.cpp.i

parameters.s: parameters.cpp.s
.PHONY : parameters.s

# target to generate assembly for a file
parameters.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/parameters.cpp.s
.PHONY : parameters.cpp.s

pdbio.o: pdbio.cpp.o
.PHONY : pdbio.o

# target to build an object file
pdbio.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/pdbio.cpp.o
.PHONY : pdbio.cpp.o

pdbio.i: pdbio.cpp.i
.PHONY : pdbio.i

# target to preprocess a source file
pdbio.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/pdbio.cpp.i
.PHONY : pdbio.cpp.i

pdbio.s: pdbio.cpp.s
.PHONY : pdbio.s

# target to generate assembly for a file
pdbio.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/pdbio.cpp.s
.PHONY : pdbio.cpp.s

timer.o: timer.cpp.o
.PHONY : timer.o

# target to build an object file
timer.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/timer.cpp.o
.PHONY : timer.cpp.o

timer.i: timer.cpp.i
.PHONY : timer.i

# target to preprocess a source file
timer.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/timer.cpp.i
.PHONY : timer.cpp.i

timer.s: timer.cpp.s
.PHONY : timer.s

# target to generate assembly for a file
timer.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/timer.cpp.s
.PHONY : timer.cpp.s

wrapper.o: wrapper.cpp.o
.PHONY : wrapper.o

# target to build an object file
wrapper.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/wrapper.cpp.o
.PHONY : wrapper.cpp.o

wrapper.i: wrapper.cpp.i
.PHONY : wrapper.i

# target to preprocess a source file
wrapper.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/wrapper.cpp.i
.PHONY : wrapper.cpp.i

wrapper.s: wrapper.cpp.s
.PHONY : wrapper.s

# target to generate assembly for a file
wrapper.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/wrapper.cpp.s
.PHONY : wrapper.cpp.s

xyzio.o: xyzio.cpp.o
.PHONY : xyzio.o

# target to build an object file
xyzio.cpp.o:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/xyzio.cpp.o
.PHONY : xyzio.cpp.o

xyzio.i: xyzio.cpp.i
.PHONY : xyzio.i

# target to preprocess a source file
xyzio.cpp.i:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/xyzio.cpp.i
.PHONY : xyzio.cpp.i

xyzio.s: xyzio.cpp.s
.PHONY : xyzio.s

# target to generate assembly for a file
xyzio.cpp.s:
	$(MAKE) -f CMakeFiles/mt.dir/build.make CMakeFiles/mt.dir/xyzio.cpp.s
.PHONY : xyzio.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... mt"
	@echo "... rebuild_cache"
	@echo "... configreader.o"
	@echo "... configreader.i"
	@echo "... configreader.s"
	@echo "... dcdio.o"
	@echo "... dcdio.i"
	@echo "... dcdio.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... num_test.o"
	@echo "... num_test.i"
	@echo "... num_test.s"
	@echo "... parameters.o"
	@echo "... parameters.i"
	@echo "... parameters.s"
	@echo "... pdbio.o"
	@echo "... pdbio.i"
	@echo "... pdbio.s"
	@echo "... timer.o"
	@echo "... timer.i"
	@echo "... timer.s"
	@echo "... wrapper.o"
	@echo "... wrapper.i"
	@echo "... wrapper.s"
	@echo "... xyzio.o"
	@echo "... xyzio.i"
	@echo "... xyzio.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

