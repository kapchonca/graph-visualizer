cmake_minimum_required(VERSION 3.5)

project(lab4)

add_executable(lab4 main.cc bmpHandler.cc graphHandler.cc graphVisualizator.cc graphMath.cc)

set(CMAKE_BUILD_TYPE Release)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -std=c++17 -g")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})