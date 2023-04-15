cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)

project(sdf LANGUAGES CXX)

add_executable(${CMAKE_PROJECT_NAME} main.cpp include/stb_image.h include/stb_image_write.h)