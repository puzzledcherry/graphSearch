cmake_minimum_required(VERSION 3.23)
project(graphsProjRevise)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(graphsProjRevise
        graph.cpp
        graph.h
        graphtest.cpp
        main.cpp)
