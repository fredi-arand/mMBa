cmake_minimum_required(VERSION 3.0)
project(mMBa)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

add_compile_options(-Wall -Wextra -Werror)

add_executable(mMBa src/main.cpp)
target_include_directories(mMBa PRIVATE src)
target_compile_features(mMBa PRIVATE cxx_std_17)

target_link_libraries(mMBa ${CONAN_LIBS})
