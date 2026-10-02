#!/bin/sh
# Requires CMake and Eigen, e.g. `brew install cmake eigen`
cmake -S . -B build && cmake --build build --config Release -j
