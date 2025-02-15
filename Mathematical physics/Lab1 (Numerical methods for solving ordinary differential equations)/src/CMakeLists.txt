cmake_minimum_required(VERSION 3.29)
project(Lab1)

set(CMAKE_CXX_STANDARD 20)

add_executable(Lab1 main.cpp
        runge_kutta_4th_order.h
        "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
        euler.h
        adams_bashforth_4th_order.h)