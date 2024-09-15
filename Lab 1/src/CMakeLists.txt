cmake_minimum_required(VERSION 3.27)
project(Lab_1)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(Lab_1
        Gauss.cpp
        Gauss.h
        main.cpp
        "QR decomposition.cpp"
        "QR decomposition.h")
