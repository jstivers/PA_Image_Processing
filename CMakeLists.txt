cmake_minimum_required(VERSION 3.6)
project(PA_Image_Processing)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES filter.cpp bmplib.cpp tommy.bmp usc_ucla_wikimedia.bmp daily-trojan.bmp bw_tile_wikimedia.bmp)
add_executable( lab_bmplib ${SOURCE_FILES})


