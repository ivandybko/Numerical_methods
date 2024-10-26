# Lab 2: Iterative methods for solving systems of linear algebraic equations

This repository implements several iterative methods for solving systems of linear equations. The methods include the Simple Iteration Method, Jacobi Method, Gauss — Seidel Method and Successive Over-Relaxation (SOR). Additionally, utilities for handling stopping criteria and solutions are provided.

## Files
- **[`SimpleIterationMethod.h`](Lab2/src/SimpleIterationMethod.h)**: Implements the Simple Iteration Method, where the system is transformed to the form ![Equation](https://latex.codecogs.com/svg.latex?x%20%3D%20Bx%20%2B%20c), and the iteration proceeds. This method can optionally use a Jacobi preconditioner.

- **[`JacobiMethod.h`](Lab2/src/JacobiMethod.h)**: Implements the Jacobi Method, an iterative method that updates all unknowns simultaneously using values from the previous iteration.
  
- **[`SeidelMethod.h`](Lab2/src/SeidelMethod.h)**: Implements the Gauss-Seidel Method, an iterative method using the latest computed values to speed up convergence.
  
- **[`SuccessiveRelaxationMethod.h`](Lab2/src/SuccessiveRelaxationMethod.h)**: Contains the implementation of the Successive Over-Relaxation (SOR) Method, which accelerates the Gauss-Seidel method by introducing a relaxation factor.

- **[`StoppingCriteria.h`](Lab2/src/StoppingCriteria.h)**: Defines stopping criteria to control the iteration process based on a given tolerance.

- **[`Solution.h`](Lab2/src/Solution.h)**: A utility file that defines structures and functions for managing solutions to the system of equations.

The Gauss-Seidel and Successive Over-Relaxation (SOR) methods have been adapted for solving systems with tridiagonal matrices, which means that they are optimized for matrices where non-zero elements only appear on the main diagonal, the diagonal directly above, and the diagonal directly below. This structure allows for more efficient memory usage and computation by taking advantage of the tridiagonal form, leading to faster convergence.

## Example Usage

The following is an example of how to use the implemented methods to solve systems of linear equations:

```cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include "SimpleIterationMethod.h"
#include "JacobiMethod.h"
#include "SeidelMethod.h"
#include "SuccessiveRelaxationMethod.h"

int main() {
    std::cout << "Enter the directory path: ";
    std::string path;
    std::cin >> path;
    
    // Reading data from files
    auto data = *readData<double>(path);

    // Solving the system using different methods
    auto sol_simpleiter = SimpleIterationMethod<double>(data, 100000, 1e-2, 1e-7);
    if (!sol_simpleiter.x.empty()) {
        std::cout << "Solution found in " << sol_simpleiter.iter 
                  << " iterations using Simple Iteration Method (C norm = " 
                  << sol_simpleiter.norm_C << "): " << sol_simpleiter.x << '\n';
    }

    auto sol_jacobi = JacobiMethod<double>(data, 100000, 1e-2);
    if (!sol_jacobi.x.empty()) {
        std::cout << "Solution found in " << sol_jacobi.iter 
                  << " iterations using Jacobi Method (C norm = " 
                  << sol_jacobi.norm_C << "): " << sol_jacobi.x << '\n';
    }

    auto sol_seidel = SeidelMethod<double>(data, 100000, 1e-2);
    if (!sol_seidel.x.empty()) {
        std::cout << "Solution found in " << sol_seidel.iter 
                  << " iterations using Gauss-Seidel Method (C norm = " 
                  << sol_seidel.norm_C << "): " << sol_seidel.x << '\n';
    }

    auto sol_relaxation = SuccessiveRelaxationMethod<double>(data, 10000, 1e-2, 0.4);
    if (!sol_relaxation.x.empty()) {
        std::cout << "Solution found in " << sol_relaxation.iter 
                  << " iterations using SOR Method (C norm = " 
                  << sol_relaxation.norm_C << "): " << sol_relaxation.x << '\n';
    }

    auto sol_relaxation_tridiag = SuccessiveRelaxationMethod<double>(
                                                                      *std::get<0>(tridiagonal_matrix),
                                                                      *std::get<1>(tridiagonal_matrix),
                                                                      *std::get<2>(tridiagonal_matrix),
                                                                      *std::get<3>(tridiagonal_matrix),
                                                                      10000, 1e-4, 1.5);
    if (!sol_relaxation_tridiag.x.empty()) {
        std::cout << "Solution found in " << sol_relaxation_tridiag.iter << " iterations using the Successive Over-Relaxation method for a tridiagonal matrix: " << sol_relaxation_tridiag.x << '\n';
    }
}
```
## Dependencies

To build and run this project, ensure you have the following dependencies installed:

- **C++11 or later**: The code uses modern C++ features, so you need a compatible compiler (e.g., GCC, Clang, or MSVC).
  
- **Matrix Operations (from Lab 1)**: The project relies on matrix operations implemented in the file `matrix_operations.h` from the previous lab. Make sure this file is included in your project’s directory structure at the specified path.
