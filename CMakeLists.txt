cmake_minimum_required(VERSION 3.16)
project(mMBa CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)

add_executable(mMBa
  src/main.cpp
  src/CommandLine.cpp
  src/DistanceField.cpp
  src/PoreMorphology.cpp
  src/VoxelVolume.cpp
)
target_compile_features(mMBa PRIVATE cxx_std_17)
target_link_libraries(mMBa PRIVATE Eigen3::Eigen Threads::Threads)
