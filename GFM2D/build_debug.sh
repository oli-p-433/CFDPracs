#!/bin/bash

# Variables
COMPILER=g++
DEBUG_FLAGS="-g -O0"          # -g for debug info, -O0 to disable optimizations
OUTPUT_EXECUTABLE="fvm" # Name of the output executable
SOURCE_FILES="solver.C fvm.C"  # Source files
INCLUDE_PATHS="-I./include"    # If you have any include directories, add them here
LIBRARIES=""                   # If you need to link any libraries, specify them here

# Compile and link
$COMPILER $DEBUG_FLAGS $INCLUDE_PATHS $SOURCE_FILES $LIBRARIES -o $OUTPUT_EXECUTABLE -Wall -Wextra -pedantic -fsanitize=address

# Check if the build was successful
if [ $? -eq 0 ]; then
    echo "Build successful. Run the following command to debug:"
    echo "gdb ./$OUTPUT_EXECUTABLE"
else
    echo "Build failed. Please check for errors."
fi
