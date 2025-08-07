# This file is part of FrontUQ-2024-proceedings-likelihood-ratios.
# Copyright (C) 2025 Emil Loevbak emil.loevbak@kit.edu
#
# FrontUQ-2024-proceedings-likelihood-ratios is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FrontUQ-2024-proceedings-likelihood-ratios is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FrontUQ-2024-proceedings-likelihood-ratios.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.14)
project(FrontUQ2024LikelihoodRatios LANGUAGES CXX)

# Set compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -O3 -DNDEBUG -fopenmp")

# Add your executable (replace main.cpp with your source files)
add_executable(frontuq_likelihood_ratios likelihoodRatioMC.cpp heat_equation_mc.cpp heat_equation_grid.cpp PeriodicTridiagonalMatrix.cpp)