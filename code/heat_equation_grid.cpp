/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#include "heat_equation_grid.hpp"
#include <iostream>
#include <cassert>

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

        std::cout << "Solving heat equation with diffusion coefficient: " << diffusionCoefficient << std::endl;

        // Do time stepping
        for (size_t k = 0; k * dt < endTime; ++k)
        {
            // Solve system
            matrix.solve(*solution, *solution_tmp);
            std::cout << "Time step " << k << std::endl;
            std::cout << "Solution_tmp: ";
            for (const auto &elem : *solution_tmp)
            {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
            std::cout << "Solution: ";
            for (const auto &elem : *solution)
            {
                std::cout << elem << " ";
            }
            std::cout << std::endl;
            std::swap(solution, solution_tmp);
        }
        assert(false);
        return *solution;
    }

}; // namespace solvers