/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#pragma once
#include <vector>

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
        std::vector<double> solution;

    public:
        HeatEquationMonteCarlo(double domainLength, size_t numberOfCells, double dt, double endTime, size_t numberOfParticles);
        std::vector<double> solve(double diffusionCoefficient);
    };
}; // namespace solvers