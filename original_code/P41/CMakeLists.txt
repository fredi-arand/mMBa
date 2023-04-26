project(P41)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
#add_executable(${PROJECT_NAME} ${SRC_LIST})

#set(CMAKE_C_COMPILER "/usr/local/bin/gcc-8")
#set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-8")

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ../eigen-new)
add_definitions( -DNOMINMAX )
add_executable(morphology_statistics main.cpp voxelvolume.h voxelio.h definitions.h
  morphologystatistics.h morphologystatistics.cpp logarithmic_bins.h histogram.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -fopenmp -fconcepts")
set(CMAKE_BUILD_TYPE release)
