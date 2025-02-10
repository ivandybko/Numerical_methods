cmake_minimum_required(VERSION 3.29)
project(src)

set(CMAKE_CXX_STANDARD 20)

add_executable(src main.cpp
        FindRoots.h
        "../../Lab3 (Solving interpolation problems)/src/grid_generators.h"
        "../../Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h")
