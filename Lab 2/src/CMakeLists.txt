cmake_minimum_required(VERSION 3.27)
project(Lab_2)

set(CMAKE_CXX_STANDARD 17)

add_executable(Lab_2 main.cpp
        "JacobiMethod.h"
        "SimpleIterationMethod.h"
        SeidelMethod.h
        "../../Lab 1/src/matrix_operations.h")
