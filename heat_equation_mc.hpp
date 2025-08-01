/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#pragma once

#include <Eigen/Dense>

namespace solvers
{
    class HeatEquationMonteCarlo
    {
    private:
        double domainLength;
        size_t numberOfCells;
        double dt;
        double endTime;
        size_t numberOfParticles;
        Eigen::VectorXd solution;

    public:
        HeatEquationMonteCarlo(double domainLength, size_t numberOfCells, double dt, double endTime, size_t numberOfParticles);
        Eigen::VectorXd solve(double diffusionCoefficient);
    };
}; // namespace solvers