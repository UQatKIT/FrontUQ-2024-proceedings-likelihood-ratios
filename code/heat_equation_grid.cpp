/*
 * This file is part of FrontUQ-2024-proceedings-likelihood-ratios.
 * Copyright (C) 2025 Emil Loevbak emil.loevbak@kit.edu
 *
 * FrontUQ-2024-proceedings-likelihood-ratios is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FrontUQ-2024-proceedings-likelihood-ratios is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FrontUQ-2024-proceedings-likelihood-ratios.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "heat_equation_grid.hpp"
#include <iostream>

namespace solvers
{

    /**
     * Constructor
     */
    HeatEquationGrid::HeatEquationGrid(double domainLength, size_t numberOfCells, double dt, double endTime)
        : matrix(numberOfCells), domainLength(domainLength), numberOfCells(numberOfCells), dt(dt), endTime(endTime)
    {
        // Initialize solution vectors
        solution = std::make_shared<std::vector<double>>(numberOfCells, 0.0);
        solution_tmp = std::make_shared<std::vector<double>>(numberOfCells, 0.0);
    }

    std::vector<double> HeatEquationGrid::solve(double diffusionCoefficient)
    {
        // Set up matrix system
        double dx = domainLength / numberOfCells;
        double alpha = diffusionCoefficient * dt / (dx * dx);

        // Main diagonal
        for (size_t i = 0; i < numberOfCells; i++)
        {
            matrix(i, i) = 1 + 2 * alpha;
        }
        // Off-diagonals
        for (size_t i = 0; i < numberOfCells - 1; i++)
        {
            matrix(i, i + 1) = -alpha;
            matrix(i + 1, i) = -alpha;
        }
        // Periodic boundaries
        matrix(0, numberOfCells - 1) = -alpha;
        matrix(numberOfCells - 1, 0) = -alpha;

        // Initial condition = dirac delta at left boundary
        // Solution vector represents density in the given cell
        std::fill(solution->begin(), solution->end(), 0.0);
        (*solution)[0] = 1.0 / dx;

        // Do time stepping
        for (size_t k = 0; k * dt < endTime; ++k)
        {
            // Solve system
            matrix.solve(*solution, *solution_tmp);
            std::swap(solution, solution_tmp);
        }
        return *solution;
    }

}; // namespace solvers