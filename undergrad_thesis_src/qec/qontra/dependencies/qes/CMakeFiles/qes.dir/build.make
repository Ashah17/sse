# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.4/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/aryan/Desktop/quantumResearch/qec/qontra

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/aryan/Desktop/quantumResearch/qec/qontra

# Include any dependencies generated for this target.
include dependencies/qes/CMakeFiles/qes.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dependencies/qes/CMakeFiles/qes.dir/compiler_depend.make

# Include the progress variables for this target.
include dependencies/qes/CMakeFiles/qes.dir/progress.make

# Include the compile flags for this target's objects.
include dependencies/qes/CMakeFiles/qes.dir/flags.make

dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o: dependencies/qes/CMakeFiles/qes.dir/flags.make
dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o: dependencies/qes/src/qes/lang/parse.cpp
dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o: dependencies/qes/CMakeFiles/qes.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/aryan/Desktop/quantumResearch/qec/qontra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o -MF CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o.d -o CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o -c /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/lang/parse.cpp

dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/qes.dir/src/qes/lang/parse.cpp.i"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/lang/parse.cpp > CMakeFiles/qes.dir/src/qes/lang/parse.cpp.i

dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/qes.dir/src/qes/lang/parse.cpp.s"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/lang/parse.cpp -o CMakeFiles/qes.dir/src/qes/lang/parse.cpp.s

dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o: dependencies/qes/CMakeFiles/qes.dir/flags.make
dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o: dependencies/qes/src/qes/lang/parse_impl.cpp
dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o: dependencies/qes/CMakeFiles/qes.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/aryan/Desktop/quantumResearch/qec/qontra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o -MF CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o.d -o CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o -c /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/lang/parse_impl.cpp

dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.i"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/lang/parse_impl.cpp > CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.i

dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.s"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/lang/parse_impl.cpp -o CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.s

dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o: dependencies/qes/CMakeFiles/qes.dir/flags.make
dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o: dependencies/qes/src/qes/util/lexer.cpp
dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o: dependencies/qes/CMakeFiles/qes.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/aryan/Desktop/quantumResearch/qec/qontra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o -MF CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o.d -o CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o -c /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/util/lexer.cpp

dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/qes.dir/src/qes/util/lexer.cpp.i"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/util/lexer.cpp > CMakeFiles/qes.dir/src/qes/util/lexer.cpp.i

dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/qes.dir/src/qes/util/lexer.cpp.s"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/util/lexer.cpp -o CMakeFiles/qes.dir/src/qes/util/lexer.cpp.s

dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o: dependencies/qes/CMakeFiles/qes.dir/flags.make
dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o: dependencies/qes/src/qes/util/llparser.cpp
dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o: dependencies/qes/CMakeFiles/qes.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/aryan/Desktop/quantumResearch/qec/qontra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o -MF CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o.d -o CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o -c /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/util/llparser.cpp

dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/qes.dir/src/qes/util/llparser.cpp.i"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/util/llparser.cpp > CMakeFiles/qes.dir/src/qes/util/llparser.cpp.i

dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/qes.dir/src/qes/util/llparser.cpp.s"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && /opt/homebrew/bin/g++-13 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/src/qes/util/llparser.cpp -o CMakeFiles/qes.dir/src/qes/util/llparser.cpp.s

# Object files for target qes
qes_OBJECTS = \
"CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o" \
"CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o" \
"CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o" \
"CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o"

# External object files for target qes
qes_EXTERNAL_OBJECTS =

dependencies/qes/libqes.a: dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse.cpp.o
dependencies/qes/libqes.a: dependencies/qes/CMakeFiles/qes.dir/src/qes/lang/parse_impl.cpp.o
dependencies/qes/libqes.a: dependencies/qes/CMakeFiles/qes.dir/src/qes/util/lexer.cpp.o
dependencies/qes/libqes.a: dependencies/qes/CMakeFiles/qes.dir/src/qes/util/llparser.cpp.o
dependencies/qes/libqes.a: dependencies/qes/CMakeFiles/qes.dir/build.make
dependencies/qes/libqes.a: dependencies/qes/CMakeFiles/qes.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/aryan/Desktop/quantumResearch/qec/qontra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libqes.a"
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && $(CMAKE_COMMAND) -P CMakeFiles/qes.dir/cmake_clean_target.cmake
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/qes.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dependencies/qes/CMakeFiles/qes.dir/build: dependencies/qes/libqes.a
.PHONY : dependencies/qes/CMakeFiles/qes.dir/build

dependencies/qes/CMakeFiles/qes.dir/clean:
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes && $(CMAKE_COMMAND) -P CMakeFiles/qes.dir/cmake_clean.cmake
.PHONY : dependencies/qes/CMakeFiles/qes.dir/clean

dependencies/qes/CMakeFiles/qes.dir/depend:
	cd /Users/aryan/Desktop/quantumResearch/qec/qontra && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/aryan/Desktop/quantumResearch/qec/qontra /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes /Users/aryan/Desktop/quantumResearch/qec/qontra /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes /Users/aryan/Desktop/quantumResearch/qec/qontra/dependencies/qes/CMakeFiles/qes.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dependencies/qes/CMakeFiles/qes.dir/depend

