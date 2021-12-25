project(P32)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
#add_executable(${PROJECT_NAME} ${SRC_LIST})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ../eigen-new)
add_definitions( -DNOMINMAX )
add_executable(pore_morphology main.cpp voxelvolume.h distancefield.h
  poremorphology.h examplespaper.h importhomberg.h thesis_helpers.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fopenmp -O3")
set(CMAKE_BUILD_TYPE release)
