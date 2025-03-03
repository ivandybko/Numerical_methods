cmake_minimum_required(VERSION 3.30)
project(src)

set(CMAKE_CXX_STANDARD 20)

add_executable(src main.cpp
        heat_1d_eq_spatial_k.h)
