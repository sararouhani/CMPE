cmake_minimum_required(VERSION 3.10)
project(CMPE)
set(CMAKE_CXX_STANDARD 11)
set(default_build_type "Release")

# Include ortools lib and include directories here
#include_directories(software/ortools/include)
#link_directories(software/ortools/lib)

add_definitions(-DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(CMPE main.cpp MN.cpp Knapsack.cpp)
add_executable(generate_mps main_generate_mps.cpp MN.cpp Knapsack.cpp)
add_executable(generate_q main_generate_q.cpp MN.cpp Knapsack.cpp)
add_executable(MCKP_Greedy main_mckp_test.cpp Knapsack.cpp)

target_link_libraries(generate_mps protobuf glog gflags CbcSolver Cbc OsiCbc Cgl ClpSolver Clp OsiClp Osi CoinUtils ortools)
target_link_libraries(MCKP_Greedy protobuf glog gflags CbcSolver Cbc OsiCbc Cgl ClpSolver Clp OsiClp Osi CoinUtils ortools)
