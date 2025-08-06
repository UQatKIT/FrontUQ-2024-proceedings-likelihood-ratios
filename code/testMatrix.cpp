#include "PeriodicTridiagonalMatrix.hpp"
#include <iostream>
#include <vector>

using namespace solvers;

int main()
{
    // Matrix size
    const int n = 5;

    // Build the 3x3 matrix: 2 on the diagonal, 1 elsewhere
    PeriodicTridiagonalMatrix mat(n);

    mat(0, 0) = 2.0;
    mat(0, 1) = 1.0;
    mat(0, 4) = 1.0;
    mat(1, 0) = 1.0;
    mat(1, 1) = 2.0;
    mat(1, 2) = 1.0;
    mat(2, 1) = 2.0;
    mat(2, 2) = 4.0;
    mat(2, 3) = 1.0;
    mat(3, 2) = 1.0;
    mat(3, 3) = 3.0;
    mat(3, 4) = 2.0;
    mat(4, 0) = 1.0;
    mat(4, 3) = 2.0;
    mat(4, 4) = 5.0;

    // Right-hand side vector
    std::vector<double> rhs = {1.0, 2.0, 3.0, 4.0, 5.0};

    // Solution vector
    std::vector<double> x(n);

    // Solve the system
    mat.solve(rhs, x);

    // Print the solution
    std::cout << "Solution:" << std::endl;
    for (double xi : x)
    {
        std::cout << xi << " ";
    }
    std::cout << std::endl;

    return 0;
}