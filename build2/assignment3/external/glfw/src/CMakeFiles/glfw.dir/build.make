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
include assignment3/external/glfw/src/CMakeFiles/glfw.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.make

# Include the progress variables for this target.
include assignment3/external/glfw/src/CMakeFiles/glfw.dir/progress.make

# Include the compile flags for this target's objects.
include assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make

assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.o: ../assignment3/external/glfw/src/context.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.o -MF CMakeFiles/glfw.dir/context.c.o.d -o CMakeFiles/glfw.dir/context.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/context.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/context.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/context.c > CMakeFiles/glfw.dir/context.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/context.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/context.c -o CMakeFiles/glfw.dir/context.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.o: ../assignment3/external/glfw/src/init.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.o -MF CMakeFiles/glfw.dir/init.c.o.d -o CMakeFiles/glfw.dir/init.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/init.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/init.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/init.c > CMakeFiles/glfw.dir/init.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/init.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/init.c -o CMakeFiles/glfw.dir/init.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.o: ../assignment3/external/glfw/src/input.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.o -MF CMakeFiles/glfw.dir/input.c.o.d -o CMakeFiles/glfw.dir/input.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/input.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/input.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/input.c > CMakeFiles/glfw.dir/input.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/input.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/input.c -o CMakeFiles/glfw.dir/input.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.o: ../assignment3/external/glfw/src/monitor.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.o -MF CMakeFiles/glfw.dir/monitor.c.o.d -o CMakeFiles/glfw.dir/monitor.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/monitor.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/monitor.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/monitor.c > CMakeFiles/glfw.dir/monitor.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/monitor.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/monitor.c -o CMakeFiles/glfw.dir/monitor.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o: ../assignment3/external/glfw/src/vulkan.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o -MF CMakeFiles/glfw.dir/vulkan.c.o.d -o CMakeFiles/glfw.dir/vulkan.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/vulkan.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/vulkan.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/vulkan.c > CMakeFiles/glfw.dir/vulkan.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/vulkan.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/vulkan.c -o CMakeFiles/glfw.dir/vulkan.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.o: ../assignment3/external/glfw/src/window.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.o -MF CMakeFiles/glfw.dir/window.c.o.d -o CMakeFiles/glfw.dir/window.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/window.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/window.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/window.c > CMakeFiles/glfw.dir/window.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/window.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -Wdeclaration-after-statement -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/window.c -o CMakeFiles/glfw.dir/window.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o: ../assignment3/external/glfw/src/cocoa_init.m
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o -MF CMakeFiles/glfw.dir/cocoa_init.m.o.d -o CMakeFiles/glfw.dir/cocoa_init.m.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_init.m

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_init.m.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_init.m > CMakeFiles/glfw.dir/cocoa_init.m.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_init.m.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_init.m -o CMakeFiles/glfw.dir/cocoa_init.m.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o: ../assignment3/external/glfw/src/cocoa_joystick.m
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o -MF CMakeFiles/glfw.dir/cocoa_joystick.m.o.d -o CMakeFiles/glfw.dir/cocoa_joystick.m.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_joystick.m

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_joystick.m.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_joystick.m > CMakeFiles/glfw.dir/cocoa_joystick.m.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_joystick.m.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_joystick.m -o CMakeFiles/glfw.dir/cocoa_joystick.m.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o: ../assignment3/external/glfw/src/cocoa_monitor.m
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o -MF CMakeFiles/glfw.dir/cocoa_monitor.m.o.d -o CMakeFiles/glfw.dir/cocoa_monitor.m.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_monitor.m

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_monitor.m.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_monitor.m > CMakeFiles/glfw.dir/cocoa_monitor.m.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_monitor.m.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_monitor.m -o CMakeFiles/glfw.dir/cocoa_monitor.m.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o: ../assignment3/external/glfw/src/cocoa_window.m
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o -MF CMakeFiles/glfw.dir/cocoa_window.m.o.d -o CMakeFiles/glfw.dir/cocoa_window.m.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_window.m

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_window.m.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_window.m > CMakeFiles/glfw.dir/cocoa_window.m.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_window.m.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_window.m -o CMakeFiles/glfw.dir/cocoa_window.m.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o: ../assignment3/external/glfw/src/cocoa_time.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o -MF CMakeFiles/glfw.dir/cocoa_time.c.o.d -o CMakeFiles/glfw.dir/cocoa_time.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_time.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/cocoa_time.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_time.c > CMakeFiles/glfw.dir/cocoa_time.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/cocoa_time.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/cocoa_time.c -o CMakeFiles/glfw.dir/cocoa_time.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o: ../assignment3/external/glfw/src/posix_thread.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o -MF CMakeFiles/glfw.dir/posix_thread.c.o.d -o CMakeFiles/glfw.dir/posix_thread.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/posix_thread.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/posix_thread.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/posix_thread.c > CMakeFiles/glfw.dir/posix_thread.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/posix_thread.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/posix_thread.c -o CMakeFiles/glfw.dir/posix_thread.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o: ../assignment3/external/glfw/src/nsgl_context.m
assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o -MF CMakeFiles/glfw.dir/nsgl_context.m.o.d -o CMakeFiles/glfw.dir/nsgl_context.m.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/nsgl_context.m

assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/nsgl_context.m.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/nsgl_context.m > CMakeFiles/glfw.dir/nsgl_context.m.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/nsgl_context.m.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/nsgl_context.m -o CMakeFiles/glfw.dir/nsgl_context.m.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o: ../assignment3/external/glfw/src/egl_context.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o -MF CMakeFiles/glfw.dir/egl_context.c.o.d -o CMakeFiles/glfw.dir/egl_context.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/egl_context.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/egl_context.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/egl_context.c > CMakeFiles/glfw.dir/egl_context.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/egl_context.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/egl_context.c -o CMakeFiles/glfw.dir/egl_context.c.s

assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/flags.make
assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o: ../assignment3/external/glfw/src/osmesa_context.c
assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o: assignment3/external/glfw/src/CMakeFiles/glfw.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o -MF CMakeFiles/glfw.dir/osmesa_context.c.o.d -o CMakeFiles/glfw.dir/osmesa_context.c.o -c /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/osmesa_context.c

assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glfw.dir/osmesa_context.c.i"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/osmesa_context.c > CMakeFiles/glfw.dir/osmesa_context.c.i

assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glfw.dir/osmesa_context.c.s"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src/osmesa_context.c -o CMakeFiles/glfw.dir/osmesa_context.c.s

# Object files for target glfw
glfw_OBJECTS = \
"CMakeFiles/glfw.dir/context.c.o" \
"CMakeFiles/glfw.dir/init.c.o" \
"CMakeFiles/glfw.dir/input.c.o" \
"CMakeFiles/glfw.dir/monitor.c.o" \
"CMakeFiles/glfw.dir/vulkan.c.o" \
"CMakeFiles/glfw.dir/window.c.o" \
"CMakeFiles/glfw.dir/cocoa_init.m.o" \
"CMakeFiles/glfw.dir/cocoa_joystick.m.o" \
"CMakeFiles/glfw.dir/cocoa_monitor.m.o" \
"CMakeFiles/glfw.dir/cocoa_window.m.o" \
"CMakeFiles/glfw.dir/cocoa_time.c.o" \
"CMakeFiles/glfw.dir/posix_thread.c.o" \
"CMakeFiles/glfw.dir/nsgl_context.m.o" \
"CMakeFiles/glfw.dir/egl_context.c.o" \
"CMakeFiles/glfw.dir/osmesa_context.c.o"

# External object files for target glfw
glfw_EXTERNAL_OBJECTS =

assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/context.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/init.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/input.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/monitor.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/vulkan.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/window.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_init.m.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_joystick.m.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_monitor.m.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_window.m.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/cocoa_time.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/posix_thread.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/nsgl_context.m.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/egl_context.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/osmesa_context.c.o
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/build.make
assignment3/external/glfw/src/libglfw3.a: assignment3/external/glfw/src/CMakeFiles/glfw.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/apple/documents/CompFabAssignment/build2/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Linking C static library libglfw3.a"
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw.dir/cmake_clean_target.cmake
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glfw.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
assignment3/external/glfw/src/CMakeFiles/glfw.dir/build: assignment3/external/glfw/src/libglfw3.a
.PHONY : assignment3/external/glfw/src/CMakeFiles/glfw.dir/build

assignment3/external/glfw/src/CMakeFiles/glfw.dir/clean:
	cd /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src && $(CMAKE_COMMAND) -P CMakeFiles/glfw.dir/cmake_clean.cmake
.PHONY : assignment3/external/glfw/src/CMakeFiles/glfw.dir/clean

assignment3/external/glfw/src/CMakeFiles/glfw.dir/depend:
	cd /Users/apple/documents/CompFabAssignment/build2 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/apple/documents/CompFabAssignment /Users/apple/documents/CompFabAssignment/assignment3/external/glfw/src /Users/apple/documents/CompFabAssignment/build2 /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src /Users/apple/documents/CompFabAssignment/build2/assignment3/external/glfw/src/CMakeFiles/glfw.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : assignment3/external/glfw/src/CMakeFiles/glfw.dir/depend
