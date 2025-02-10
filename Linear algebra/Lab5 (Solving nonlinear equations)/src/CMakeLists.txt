cmake_minimum_required(VERSION 3.29)
project(src)

set(CMAKE_CXX_STANDARD 20)

add_executable(src main.cpp
        FindRoots.h
        ../../Lab3/src/grid_generators.h
        ../../Lab1/src/matrix_operations.h)
