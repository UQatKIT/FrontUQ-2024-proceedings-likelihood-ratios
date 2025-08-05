/*
 * This file is part of TODO.
 * Copyright (C) 2021 Emil Loevbak emil.loevbak@kuleuven.be
 *
 * TODO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TODO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TODO.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "TridiagonalMatrix.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

namespace solvers
{

    TridiagonalMatrix::TridiagonalMatrix(size_t dimension) : dimension(dimension)
    {
        assert(dimension > 0);

        mainDiagonal = std::make_unique<std::vector<double>>(dimension, 0.0);
        subDiagonal = std::make_unique<std::vector<double>>(dimension - 1, 0.0);
        superDiagonal = std::make_unique<std::vector<double>>(dimension - 1, 0.0);
    }

    double &TridiagonalMatrix::operator()(size_t row, size_t column)
    {
        assert(row < dimension);
        assert(column < dimension);
        int diagonal = column - row;
        assert(diagonal <= 1);

        return diagonalElement(diagonal, std::min(row, column));
    }

    double const &TridiagonalMatrix::operator()(size_t row, size_t column) const
    {
        assert(row < dimension);
        assert(column < dimension);

        int diagonal = column - row;
        if (std::abs(diagonal) <= 1)
        {
            return diagonalElement(diagonal, std::min(row, column));
        }
        else
        {
            return zero;
        }
    }

    size_t TridiagonalMatrix::size()
    {
        return dimension;
    }

    double &TridiagonalMatrix::diagonalElement(int diagonal, size_t element)
    {
        assert(std::abs(diagonal) <= 1);
        assert(element < dimension - std::abs(diagonal));

        switch (diagonal)
        {
        case -1:
            return (*subDiagonal)[element];
        case 0:
            return (*mainDiagonal)[element];
        case 1:
            return (*superDiagonal)[element];
        default:
            std::cerr << "Invalid attempt to acces a zero-defined element of the matrix as non-const!" << std::endl;
            return const_cast<double &>(zero); // This line is never executed, but the compiler complains if it is not included.
        }
    }

    // TODO: This causes a lot of branching!
    double const &TridiagonalMatrix::diagonalElement(int diagonal, size_t element) const
    {
        assert(element + std::abs(diagonal) < dimension);

        switch (diagonal)
        {
        case -1:
            return (*subDiagonal)[element];
        case 0:
            return (*mainDiagonal)[element];
        case 1:
            return (*superDiagonal)[element];
        default:
            return zero;
        }
    }

    // TODO: Can we guarantee that the matrix is diagonal dominant?
    //       If not we need to test and fall back on DGTSV.
    void TridiagonalMatrix::solve(const std::vector<double> &rhs, std::vector<double> &solution)
    {
        assert(rhs.size() == dimension);
        assert(solution.size() == dimension);

        std::vector<double> c(dimension - 1);
        std::vector<double> d(dimension - 1);

        // Implementation of the Thomas algorithm (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).
        // Based on Mathias' implementation. Note that the indices for a (the *subdiagonal here) start from 2 on Wikipedia.
        c[0] = (*superDiagonal)[0] / (*mainDiagonal)[0];
        d[0] = rhs[0] / (*mainDiagonal)[0];
        for (size_t i = 1; i < dimension - 1; ++i)
        {
            double denominator = (*mainDiagonal)[i] - (*subDiagonal)[i - 1] * c[i - 1];
            c[i] = (*superDiagonal)[i] / denominator;
            d[i] = (rhs[i] - (*subDiagonal)[i - 1] * d[i - 1]) / denominator;
        }

        solution[dimension - 1] = (rhs[dimension - 1] - (*subDiagonal)[dimension - 2] * d[dimension - 2]) /
                                  ((*mainDiagonal)[dimension - 1] - (*subDiagonal)[dimension - 2] * c[dimension - 2]);
        for (size_t i = dimension - 2; i != static_cast<size_t>(-1); --i)
        {
            solution[i] = d[i] - c[i] * solution[i + 1];
        }
    }

    void TridiagonalMatrix::transposedSolve(const std::vector<double> &rhs, std::vector<double> &solution)
    {
        // Transpose tridiagonal matrix
        std::swap(subDiagonal, superDiagonal);
        solve(rhs, solution);
        // Transpose back
        std::swap(subDiagonal, superDiagonal);
    }

    void TridiagonalMatrix::transposedMatrixVectorProductSubtract(const std::vector<double> &vector, std::vector<double> &result)
    {
        assert(result.size() == dimension);
        assert(vector.size() == dimension);
        assert(dimension > 2);

        result[0] -= (*mainDiagonal)[0] * vector[0] + (*subDiagonal)[0] * vector[1];
        for (size_t i = 1; i < dimension - 1; ++i)
        {
            result[i] -= (*superDiagonal)[i - 1] * vector[i - 1] + (*mainDiagonal)[i] * vector[i] + (*subDiagonal)[i] * vector[i + 1];
        }
        result[dimension - 1] -= (*superDiagonal)[dimension - 2] * vector[dimension - 2] + (*mainDiagonal)[dimension - 1] * vector[dimension - 1];
    }

}