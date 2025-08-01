cmake_minimum_required(VERSION 3.14)
project(FrontUQ2024LikelihoodRatios LANGUAGES CXX)

# Set compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -O3 -DNDEBUG -fopenmp")

# Download Eigen using FetchContent
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
)
FetchContent_MakeAvailable(eigen)

# Add your executable (replace main.cpp with your source files)
add_executable(frontuq_likelihood_ratios likelihoodRatioMC.cpp heat_equation_mc.cpp heat_equation_grid.cpp)

# Link Eigen
target_link_libraries(frontuq_likelihood_ratios PRIVATE Eigen3::Eigen)
target_include_directories(frontuq_likelihood_ratios PRIVATE ${eigen_SOURCE_DIR})