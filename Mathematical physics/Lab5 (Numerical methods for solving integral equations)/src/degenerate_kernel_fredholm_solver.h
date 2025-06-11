#ifndef DEGENERATE_KERNEL_FREDHOLM_SOLVER_H
#define DEGENERATE_KERNEL_FREDHOLM_SOLVER_H
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"

#include <vector>
#include <functional>

template <typename T>
std::vector<T> degenerate_kernel_fredholm_solver(T lambda, T a,	T b, std::function<T(T)> f,	std::vector<std::function<T(T)>> phis, std::vector<std::function<T(T)>> psis, int grid_size) {
	const int m = phis.size();
	T h = (b - a) / (grid_size - 1);
	std::vector<T> weights(grid_size, h);
	weights[0] /= 2;
	weights[grid_size - 1] /= 2;
	std::vector<std::vector<T>> alpha(m, std::vector<T>(m, 0.0));
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			T sum = 0.0;
			for (int k = 0; k < grid_size; ++k) {
				sum += psis[i](a + k * h) * phis[j](a + k * h) * weights[k];
			}
			alpha[i][j] = sum;
		}
	}
	std::vector<T> beta(m, 0.0);
	for (int i = 0; i < m; ++i) {
		T sum = 0.0;
		for (int k = 0; k < grid_size; ++k) {
			sum += psis[i](a + k * h) * f(a + k * h) * weights[k];
		}
		beta[i] = sum;
	}
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			if (i == j)
			{
				alpha[i][j] = 1 - lambda * alpha[i][j];
			}
			else
			{
				alpha[i][j] = -lambda * alpha[i][j];
			}
		}
	}
	std::vector<T> C = Gauss<T>({alpha, beta});
	std::vector<T> sol(grid_size);
	for (int i = 0; i < grid_size; ++i) {
		T sum = 0.0;
		for (int j = 0; j < m; ++j) {
			sum += C[j] * phis[j](a + i * h);
		}
		sol[i] = f(a + i * h) + lambda * sum;
	}
	return sol;
}

#endif //DEGENERATE_KERNEL_FREDHOLM_SOLVER_H
