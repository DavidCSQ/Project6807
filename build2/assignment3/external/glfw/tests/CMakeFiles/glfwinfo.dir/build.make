# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/apple/documents/CompFabAssignment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/apple/documents/CompFabAssignment/build2

# Include any dependencies generated for this target.
include assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/compiler_depend.make

# Include the progress variables for this target.
include assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/progress.make

# Include the compile flags for this target's objects.
include assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/flags.make

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/flags.make
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.o: ../assignment3/external/glfw/tests/glfwinfo.c
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.o -MF CMakeFiles/glfwinfo.dir/glfwinfo.c.o.d -o CMakeFiles/glfwinfo.dir/glfwinfo.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/tests/glfwinfo.c

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfwinfo.dir/glfwinfo.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/tests/glfwinfo.c > CMakeFiles/glfwinfo.dir/glfwinfo.c.i

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfwinfo.dir/glfwinfo.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/tests/glfwinfo.c -o CMakeFiles/glfwinfo.dir/glfwinfo.c.s

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/flags.make
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o: ../assignment3/external/glfw/deps/getopt.c
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o -MF CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o.d -o CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/getopt.c

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfwinfo.dir/__/deps/getopt.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/getopt.c > CMakeFiles/glfwinfo.dir/__/deps/getopt.c.i

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfwinfo.dir/__/deps/getopt.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/getopt.c -o CMakeFiles/glfwinfo.dir/__/deps/getopt.c.s

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/flags.make
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o: ../assignment3/external/glfw/deps/glad_gl.c
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o -MF CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o.d -o CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/glad_gl.c

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/glad_gl.c > CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.i

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/glad_gl.c -o CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.s

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/flags.make
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o: ../assignment3/external/glfw/deps/glad_vulkan.c
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o -MF CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o.d -o CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/glad_vulkan.c

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/glad_vulkan.c > CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.i

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/deps/glad_vulkan.c -o CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.s

# Object files for target glfwinfo
glfwinfo_OBJECTS = \
"CMakeFiles/glfwinfo.dir/glfwinfo.c.o" \
"CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o" \
"CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o" \
"CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o"

# External object files for target glfwinfo
glfwinfo_EXTERNAL_OBJECTS =

assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/glfwinfo.c.o
assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/getopt.c.o
assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_gl.c.o
assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/__/deps/glad_vulkan.c.o
assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/build.make
assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/src/libglfw3.a
assignment3/external/glfw/tests/glfwinfo: assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable glfwinfo"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glfwinfo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/build: assignment3/external/glfw/tests/glfwinfo
.PHONY : assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/build

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/clean:
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests && $(CMAKE_COMMAND) -P CMakeFiles/glfwinfo.dir/cmake_clean.cmake
.PHONY : assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/clean

assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/depend:
	cd /Users/apple/documents/CompFabAssignment/build2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/apple/documents/CompFabAssignment /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/tests /Users/apple/documents/CompFabAssignment/build2 /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : assignment3/external/glfw/tests/CMakeFiles/glfwinfo.dir/depend

