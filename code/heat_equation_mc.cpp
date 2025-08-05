/**
 * A deterministic solver for the 1D heat equation.
 * @author Emil Løvbak
 */

#include <random>
#include "heat_equation_mc.hpp"
#include <cassert>
#include <omp.h>
#include <iostream>
#include <algorithm>

namespace solvers
{

    /**
     * Constructor
     */
    HeatEquationMonteCarlo::HeatEquationMonteCarlo(double domainLength, size_t numberOfCells, double dt, double endTime, size_t numberOfParticles)
        : domainLength(domainLength), numberOfCells(numberOfCells), dt(dt), endTime(endTime), numberOfParticles(numberOfParticles) {}

    std::vector<double> HeatEquationMonteCarlo::solve(double diffusionCoefficient)
    {
        solution = std::vector<double>(numberOfCells, 0.0);
        double dx = domainLength / numberOfCells;

#pragma omp declare reduction(+ : std::vector<double> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

#pragma omp parallel
        {
            int threadNum = omp_get_thread_num();
            std::default_random_engine generator(std::random_device{}() + threadNum);
            std::normal_distribution<double> normalDistribution(0.0, sqrt(2.0 * diffusionCoefficient * dt));

#pragma omp for reduction(+ : solution)
            for (size_t p = 1; p <= numberOfParticles; ++p)
            {

                double X = dx / 2.0;

                for (size_t k = 0; k * dt < endTime; ++k)
                {
                    X += normalDistribution(generator);
                    if (X < 0.0)
                    {
                        X += domainLength * (1 + static_cast<size_t>(-X / domainLength));
                    }
                    // If not needed here as if X is in range the expression below evaluates to 0
                    X -= domainLength * (static_cast<size_t>(X / domainLength));
                }

                assert(X >= 0.0 && X < domainLength);
                size_t cellNumber = X / dx;
                solution[cellNumber] += 1.0 / numberOfParticles / dx;
            }
        }
        return solution;
    }

}; // namespace solvers