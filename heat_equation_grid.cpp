/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#include "heat_equation_grid.hpp"

namespace solvers
{

    /**
     * Constructor
     */
    HeatEquationGrid::HeatEquationGrid(double domainLength, size_t numberOfCells, double dt, double endTime)
        : matrix(numberOfCells, numberOfCells), domainLength(domainLength), numberOfCells(numberOfCells), dt(dt), endTime(endTime)
    {
        // Initialize solution vectors
        solution = std::make_shared<Eigen::VectorXd>(numberOfCells);
        solution_tmp = std::make_shared<Eigen::VectorXd>(numberOfCells);
    }

    Eigen::VectorXd HeatEquationGrid::solve(double diffusionCoefficient)
    {
        // Set up matrix system
        double dx = domainLength / numberOfCells;
        double alpha = diffusionCoefficient * dt / (dx * dx);
        matrix.setZero();

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
        solution->setZero();
        (*solution)(0) = 1.0 / dx;

        // Do time stepping
        for (size_t k = 0; k * dt < endTime; ++k)
        {
            // Solve system
            *solution_tmp = matrix.llt().solve(*solution);
            std::swap(solution, solution_tmp);
        }
        return *solution;
    }

}; // namespace solvers