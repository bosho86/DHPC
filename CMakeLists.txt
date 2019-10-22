cmake_minimum_required(VERSION 3.5.2)

project(DHPC)

##### Setting up the CXX flags #####
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEFAULT} -std=c++11 -Wall -Wextra -O3 -fno-tree-vectorize -mavx2 -mfma")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEFAULT} -std=c11 -Wall -Wextra -O3 -fno-tree-vectorize -mavx2 -mfma")



##### Include header files #####
include_directories(
        include
)
include_directories(src)

##### Source files to compile for benchmark executable #####
add_executable(sssp
        src/bellmanfordalgo.cpp
        src/dijkstra.cpp
        )


