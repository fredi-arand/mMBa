#!/bin/sh

# Setup .venv directory
python3 -m venv .venv

# Activate virtual environment
source .venv/bin/activate

# (optional) upgrade pip
# pip install --upgrade pip

# Install conan
pip install conan

# Get dependencies
mkdir build && cd build && conan install ..

# Run cmake
cmake ..

# Build the program, result should be in build/bin/mMBa
cmake --build .

# Create output dirs
mkdir ../output/
mkdir ../output/stacks

# Run program
cd bin 
./mMBa
