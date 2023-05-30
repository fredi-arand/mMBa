#!/bin/sh

# Setup .venv directory
python3 -m venv .venv

# Activate virtual environment
source .venv/bin/activate

# (optional) upgrade pip
# pip install --upgrade pip

# Install conan
pip install conan==1.60.0

# Get dependencies and build them
mkdir build && cd build && conan install --build="*" ..

# Run cmake, set Release mode
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build the program, result should be in build/bin/mMBa
cmake --build .
