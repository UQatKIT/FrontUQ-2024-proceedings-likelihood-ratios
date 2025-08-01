/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#pragma once

#include <Eigen/Dense>
#include <memory>

namespace solvers
{
    class HeatEquationGrid
    {
    private:
        Eigen::MatrixXd matrix;
        std::shared_ptr<Eigen::VectorXd> solution;
        std::shared_ptr<Eigen::VectorXd> solution_tmp;
        double domainLength;
        size_t numberOfCells;
        double dt;
        double endTime;

    public:
        HeatEquationGrid(double domainLength, size_t numberOfCells, double dt, double endTime);
        Eigen::VectorXd solve(double diffusionCoefficient);
    };
}; // namespace solvers