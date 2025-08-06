cmake_minimum_required(VERSION 3.14)
project(FrontUQ2024LikelihoodRatios LANGUAGES CXX)

# Set compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -Og -fopenmp")

# Add your executable (replace main.cpp with your source files)
add_executable(frontuq_likelihood_ratios likelihoodRatioMC.cpp heat_equation_mc.cpp heat_equation_grid.cpp PeriodicTridiagonalMatrix.cpp)

add_executable(testMatrix testMatrix.cpp PeriodicTridiagonalMatrix.cpp)