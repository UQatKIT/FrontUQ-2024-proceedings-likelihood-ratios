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

#include "PeriodicTridiagonalMatrix.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

namespace solvers
{

    PeriodicTridiagonalMatrix::PeriodicTridiagonalMatrix(size_t dimension) : dimension(dimension)
    {
        assert(dimension > 0);

        mainDiagonal = std::make_unique<std::vector<double>>(dimension, 0.0);
        subDiagonal = std::make_unique<std::vector<double>>(dimension - 1, 0.0);
        superDiagonal = std::make_unique<std::vector<double>>(dimension - 1, 0.0);
        topRight = 0.0;
        bottomLeft = 0.0;
    }

    double &PeriodicTridiagonalMatrix::operator()(size_t row, size_t column)
    {
        assert(row < dimension);
        assert(column < dimension);
        int diagonal = column - row;
        assert(abs(diagonal) <= 1 || (row == 0 && column == dimension - 1) || (row == dimension - 1 && column == 0));

        if (abs(diagonal) <= 1)
        {
            return diagonalElement(diagonal, std::min(row, column));
        }
        else if (row == 0 && column == dimension - 1)
        {
            return topRight;
        }
        else if (row == dimension - 1 && column == 0)
        {
            return bottomLeft;
        }
        else
        {
            assert(false && "Invalid diagonal index");
            std::cerr << "Invalid attempt to access a zero-defined element of the matrix!" << std::endl;
            return const_cast<double &>(zero);
        }
    }

    double const &PeriodicTridiagonalMatrix::operator()(size_t row, size_t column) const
    {
        assert(row < dimension);
        assert(column < dimension);

        int diagonal = column - row;
        if (std::abs(diagonal) <= 1)
        {
            return diagonalElement(diagonal, std::min(row, column));
        }
        else if (row == 0 && column == dimension - 1)
        {
            return topRight;
        }
        else if (row == dimension - 1 && column == 0)
        {
            return bottomLeft;
        }
        else
        {
            return zero;
        }
    }

    size_t PeriodicTridiagonalMatrix::size()
    {
        return dimension;
    }

    double &PeriodicTridiagonalMatrix::diagonalElement(int diagonal, size_t element)
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
            assert(false && "Invalid diagonal index");
            std::cerr << "Invalid attempt to acces a zero-defined element of the matrix as non-const!" << std::endl;
            return const_cast<double &>(zero);
        }
    }

    // TODO: This causes a lot of branching!
    double const &PeriodicTridiagonalMatrix::diagonalElement(int diagonal, size_t element) const
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

    void PeriodicTridiagonalMatrix::thomas(const std::vector<double> &rhs, std::vector<double> &solution)
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

    // Based on the Sherman-Morrison formula for solving a linear system with a modified matrix.
    // Formula found at (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).
    void PeriodicTridiagonalMatrix::solve(const std::vector<double> &rhs, std::vector<double> &solution)
    {
        assert(rhs.size() == dimension);
        assert(solution.size() == dimension);

        // Produce modified matrix
        double gamma = std::sqrt(bottomLeft * topRight);
        (*mainDiagonal)[0] -= gamma;
        (*mainDiagonal)[dimension - 1] -= bottomLeft * topRight / gamma;

        // Produce right hand sides
        std::vector<double> u(rhs.size(), 0.0);
        std::vector<double> v(rhs.size(), 0.0);
        u[0] = gamma;
        u[dimension - 1] = bottomLeft;
        v[0] = 1.0;
        v[dimension - 1] = topRight / gamma;

        // Solve the modified systems
        std::vector<double> q(rhs.size(), 0.0);
        std::vector<double> y(rhs.size(), 0.0);
        thomas(u, q);
        thomas(rhs, y);

        // Compute the solution
        double vy = 0.0;
        double vq = 0.0;
        std::transform(v.begin(), v.end(), y.begin(), y.begin(),
                       [](double v_elem, double y_elem)
                       { return v_elem * y_elem; });
        vy = std::accumulate(y.begin(), y.end(), 0.0);
        std::transform(v.begin(), v.end(), q.begin(), q.begin(),
                       [](double v_elem, double q_elem)
                       { return v_elem * q_elem; });
        vq = std::accumulate(q.begin(), q.end(), 0.0);
        double factor = vy / (1.0 + vq);
        std::transform(q.begin(), q.end(), y.begin(), solution.begin(),
                       [factor](double q_elem, double y_elem)
                       { return y_elem - factor * q_elem; });
    }

    void PeriodicTridiagonalMatrix::transposedSolve(const std::vector<double> &rhs, std::vector<double> &solution)
    {
        // Transpose tridiagonal matrix
        std::swap(subDiagonal, superDiagonal);
        std::swap(topRight, bottomLeft);
        solve(rhs, solution);
        // Transpose back
        std::swap(subDiagonal, superDiagonal);
        std::swap(topRight, bottomLeft);
    }

}