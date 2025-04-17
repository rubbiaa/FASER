#!/bin/bash

FASER_BASE_PATH=$1

if [ -z $FASER_BASE_PATH ]; then
	echo "Use as ./generateLib.sh /full/path/to/Faser/main/dir"
	exit 1
fi

INCLUDE_PATH="/data/sw/FASERCAL/FASER/GenFit-install/include"

# Check if the base path exists
if [ ! -d "$FASER_BASE_PATH" ]; then
    echo "Error: Directory $FASER_BASE_PATH does not exist."
    exit 1
fi

# Generate ClassesDict.cxx
echo "Generating ClassesDict.cxx..."
rootcint -f $FASER_BASE_PATH/Python_io/lib/ClassesDict.cxx -c \
    -I$INCLUDE_PATH \
    $FASER_BASE_PATH/CoreUtils/*.hh $FASER_BASE_PATH/Python_io/lib/LinkDef.h
if [ $? -ne 0 ]; then
    echo "Error: rootcint failed."
    exit 1
fi

# Compile shared library
echo "Compiling shared library..."
g++ -o $FASER_BASE_PATH/Python_io/lib/ClassesDict.so -shared -fPIC \
    `root-config --cflags --libs` \
    -I$INCLUDE_PATH \
    $FASER_BASE_PATH/CoreUtils/*.cc $FASER_BASE_PATH/Python_io/lib/ClassesDict.cxx
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

echo "Library generated successfully!"
