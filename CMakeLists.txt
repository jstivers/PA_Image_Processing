cmake_minimum_required(VERSION 3.6)
project(PA_Image_Processing)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES filter.cpp bmplib.cpp)
add_executable(lab_bmplib ${SOURCE_FILES})
