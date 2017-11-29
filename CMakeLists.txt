cmake_minimum_required(VERSION 3.9)
project(hpc_tp3)

set(CMAKE_C_STANDARD 99)

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

include_directories(${MPI_C_INCLUDE_PATH})

add_executable(distance-matrix distance-matrix.c)
target_link_libraries(distance-matrix m)
