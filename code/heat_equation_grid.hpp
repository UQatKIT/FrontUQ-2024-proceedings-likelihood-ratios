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

#include "PeriodicTridiagonalMatrix.hpp"

namespace solvers
{
    class HeatEquationGrid
    {
    private:
        PeriodicTridiagonalMatrix matrix;
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