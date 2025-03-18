cmake_minimum_required(VERSION 3.30)
project(src)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")

add_executable(src main.cpp
        heat_1d_eq_spatial_k.h
        "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h")
