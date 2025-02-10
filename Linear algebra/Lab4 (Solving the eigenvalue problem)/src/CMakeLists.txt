cmake_minimum_required(VERSION 3.29)
project(Lab3)

set(CMAKE_CXX_STANDARD 20)

add_executable(src main.cpp
        FrancisKublanovskayaQR.h
        "../../Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
        "../../Lab1 (Direct methods for solving systems of linear equations)/src/QR decomposition.h"
        "../../Lab1 (Direct methods for solving systems of linear equations)/src/Gauss.h"
        HessenbergReduction.h
        InverseIterationMethod.h
        "../../Lab2 (Iterative methods for solving systems of linear equations)/src/StoppingCriteria.h")
