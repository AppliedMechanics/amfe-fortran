# define minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR) # TODO: check version requirement

# project name and supported languages
project(amfefortran LANGUAGES CXX Fortran)

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "CMake was called")
