/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#pragma once

#include "TridiagonalMatrix.hpp"

namespace solvers
{
    class HeatEquationGrid
    {
    private:
        TridiagonalMatrix matrix;
        std::shared_ptr<std::vector<double>> solution;
        std::shared_ptr<std::vector<double>> solution_tmp;
        double domainLength;
        size_t numberOfCells;
        double dt;
        double endTime;

    public:
        HeatEquationGrid(double domainLength, size_t numberOfCells, double dt, double endTime);
        std::vector<double> solve(double diffusionCoefficient);
    };
}; // namespace solvers