cmake_minimum_required(VERSION 3.29)
project(Lab3)

set(CMAKE_CXX_STANDARD 20)

add_executable(src main.cpp
        FrancisKublanovskayaQR.h
        ../../Lab1/src/matrix_operations.h
        "../../Lab1/src/QR decomposition.h"
        ../../Lab1/src/Gauss.h
        HessenbergReduction.h
        InverseIterationMethod.h
        ../../Lab2/src/StoppingCriteria.h)
