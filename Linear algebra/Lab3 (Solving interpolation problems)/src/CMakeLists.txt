cmake_minimum_required(VERSION 3.29)
project(Lab3)

set(CMAKE_CXX_STANDARD 20)

add_executable(Lab3 main.cpp
        grid_generators.h
        NewtonInterpolation.h
        SplineInterpolation.h
        "../../Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h")
