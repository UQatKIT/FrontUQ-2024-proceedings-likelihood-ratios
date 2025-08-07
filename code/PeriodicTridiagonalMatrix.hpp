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
#include <memory>

namespace solvers
{
    /**
     * Class representing a tridiagonal matrix.
     * The class also contains functionality for solving the tridiagonal matrix system.
     */
    class PeriodicTridiagonalMatrix
    {
    public:
        /**
         * Constructor. Creates an empty matrix.
         * @param dimension The number of rows/columns of the matrix.
         */
        PeriodicTridiagonalMatrix(size_t dimension);
        /**
         * Access an element of the matrix with a given row and column index.
         * The requested element must lie on the main diagonal or the subdiagonals directly above or below the main diagonal.
         * @param row The row index.
         * @param column The column index.
         * @returns A reference to the element at the given pair of indices.
         */
        double &operator()(size_t row, size_t column);
        /**
         * Access an element of the matrix with a given row and column index.
         * This function does not have the index constraints of the non-const version, but simply requires the indices to be within bounds.
         * @param row The row index.
         * @param column The column index.
         * @returns A const reference to the element at the given pair of indices.
         */
        double const &operator()(size_t row, size_t column) const;
        /**
         * Return a given element from a given diagonal. Note that this function does not access the corner elements.
         * @param diagonal The diagonal from which to select the element. 0 corresponds with the main diagonal, -1 with the sub-diagonal and 1 with the super-diagonal.
         * @param element The element on that diagonal, starting from 0.
         * @returns A reference to the requested element on the requested diagonal.
         */
        double &diagonalElement(int diagonal, size_t element);
        /**
         * Return a given element from a given diagonal. Note that this function does not access the corner elements.
         * @param diagonal The diagonal from which to select the element. 0 corresponds with the main diagonal, -n with the n-th sub-diagonal and n with the n-th super-diagonal.
         * @param element The element on that diagonal, starting from 0.
         * @returns A const reference to the requested element on the requested diagonal.
         */
        double const &diagonalElement(int diagonal, size_t element) const;
        /**
         * Solves a linear system consisting of this matrix and the given right hand side.
         * @param rhs The given right hand side. This vector must have the same size as the number of rows/columns in the matrix.
         * @param solution The location where the solution should be stored. This vector must have the same size as the number of rows/columns in the matrix.
         */
        void solve(const std::vector<double> &rhs, std::vector<double> &solution);
        /**
         * Solves a linear system consisting of the transpose of this matrix and the given right hand side.
         * @param rhs The given right hand side. This vector must have the same size as the number of rows/columns in the matrix.
         * @param solution The location where the solution should be stored. This vector must have the same size as the number of rows/columns in the matrix.
         */
        void transposedSolve(const std::vector<double> &rhs, std::vector<double> &solution);
        /**
         * Get the number of rows/columns in the matrix.
         * @returns The number of rows/columns.
         */
        size_t size();

    private:
        /**
         * The number of rows/columns in the matrix.
         */
        size_t dimension;
        /**
         * A zero value which is used in the const operator() and diagonalElement functions.
         */
        double const zero = 0.0;
        /**
         * A vector of size dimension, representing the main diagonal of the matrix.
         */
        std::unique_ptr<std::vector<double>> mainDiagonal;
        /**
         * A vector of size dimension - 1 , representing the sub diagonal of the matrix.
         */
        std::unique_ptr<std::vector<double>> subDiagonal;
        /**
         * A vector of size dimension - 1, representing the super diagonal of the matrix.
         */
        std::unique_ptr<std::vector<double>> superDiagonal;
        /**
         * The top right element of the matrix.
         */
        double topRight;
        /**
         * The bottom left element of the matrix.
         */
        double bottomLeft;
        /**
         * Implementation of the Thomas algorithm for solving a tridiagonal matrix system. Used internally by solve.
         * @param rhs The given right hand side. This vector must have the same size as the number of rows/columns in the matrix.
         * @param solution The location where the solution should be stored. This vector must have the same size as the number of rows/columns in the matrix.
         */
        void thomas(const std::vector<double> &rhs, std::vector<double> &solution);
    };

}