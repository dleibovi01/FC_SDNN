cmake_minimum_required(VERSION 3.0.0)

# set the project name
project(Shock_Solver)

# set(CMAKE_CXX_COMPILER "/opt/ohpc/pub/compiler/intel/compiler/2021.1.1/linux/bin/intel64/icpc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl -Wall -pedantic -g -O3")
set(CMAKE_CXX_STANDARD 17)
#set(SRC_FILES "Code/*.cpp")

file(GLOB SRC_FILES Code/*.cpp)


#add the executable
add_executable(Test_shock_solver ${SRC_FILES})

