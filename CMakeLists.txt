# define minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR) # TODO: check version requirement

# project name and supported languages
project(amfefortran LANGUAGES CXX Fortran)

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


message(STATUS "Project will be installed to ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH "Installation directory for libraries")

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

add_subdirectory(amfefortran)
