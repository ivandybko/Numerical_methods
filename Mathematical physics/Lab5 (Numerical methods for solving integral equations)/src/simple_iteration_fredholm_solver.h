#ifndef SIMPLE_ITERATION_FREDHOLM_SOLVER_H
#define SIMPLE_ITERATION_FREDHOLM_SOLVER_H
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"

#include <vector>
#include <functional>
template <typename T>
std::vector<T> simple_iteration_fredholm_solve(T lambda, T a, T b, std::function<T(T, T)> K, std::function<T(T)> f, int grid_size, int max_iterations, T tolerance){
	T h = (b - a) / (grid_size - 1);
	std::vector<T> sol_old(grid_size);
	for (int i = 0; i < grid_size; ++i)
    {
		sol_old[i] = f(a + i * h);
    }
	std::vector<T> sol_new(grid_size, 0.0);
	std::vector<T> weights(grid_size, h);
	weights[0] /= 2;
	weights[grid_size - 1] /= 2;
	for (int iter = 0; iter < max_iterations; ++iter) {
		for (int i = 0; i < grid_size; ++i) {
			T sum = 0.0;
			for (int j = 0; j < grid_size; ++j) {
				sum += K(a + i * h, a + j * h) * sol_old[j] * weights[j];
			}
			sol_new[i] = f(a + i * h) + lambda * sum;
		}
		auto diff = sol_new - sol_old;
		// T max_diff = *std::max_element(diff.begin(), diff.end(), std::greater<T>());
		// if (max_diff < tolerance) {
		// 	break;
		// }
		sol_old = sol_new;
	}
    return sol_new;
}

#endif //SIMPLE_ITERATION_FREDHOLM_SOLVER_H
