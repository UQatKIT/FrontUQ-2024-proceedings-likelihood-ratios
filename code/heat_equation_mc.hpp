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