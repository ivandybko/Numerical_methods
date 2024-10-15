#ifndef LAB_2__SEIDEL_METHOD_H_
#define LAB_2__SEIDEL_METHOD_H_
#include <iostream>
#include <vector>
#include <utility>
#include "../../Lab 1/src/matrix_operations.h"
#include "Solution.h"
#include "StoppingCriteria.h"

template <typename T>
T calcNormC(const std::vector<std::vector<T>>& A){
	size_t n = A.size();
	std::vector<std::vector<T>> InvLD(n, std::vector<T>(n,0));
	for (int i = 0; i < n; ++i)
	{
		InvLD[i][i]=1/A[i][i];
		for (int j = 0; j < i; ++j)
		{
			for (int k = j; k < i; ++k)
			{
				InvLD[i][j]+=A[i][k]*InvLD[k][j];
			}
			InvLD[i][j]/=-A[i][i];
		}
	}
	std::vector<std::vector<T>> U(n, std::vector<T>(n,0));
	for (int i = 0; i < n; ++i)
	{
		for (int j = i; j < n; ++j)
		{
			U[i][j]=-A[i][j];
		}
	}
	std::vector<std::vector<T>> C(n, std::vector<T>(n,0));
	multiply(InvLD,U,C);
	std::vector<std::vector<T>> C_U(n, std::vector<T>(n,0));
	std::vector<std::vector<T>> C_L(n, std::vector<T>(n,0));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			C_L[i][j] = C[i][j];
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = i+1; j < n; ++j) {
			C_U[i][j] = C[i][j];
		}
	}
	return octahedralNorm(C_L) + octahedralNorm(C_U);
}

template <typename T>
Solution<T> SeidelMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter =1000, T eps=1e-4){
	auto *A = &data.first;
	auto *b = &data.second;
	size_t n = (*A).size();
	T normC{calcNormC<T>(*A)};
	std::vector<T> x(n, 0.0);
	for (size_t iter = 0; iter < max_iter; ++iter) {
		std::vector<T> x_old{x};
		for (size_t i = 0; i < n; ++i) {
			T sum{0.0};
			for (size_t j = 0; j < i; ++j) {
				sum += (*A)[i][j] * x[j];
			}
			for (size_t j = i + 1; j < n; ++j) {
				sum += (*A)[i][j] * x_old[j];
			}
			x[i] = ((*b)[i] - sum) / (*A)[i][i];
		}
		if (checkConvergence(x, x_old, eps)) {
			return Solution<T>{x, iter , normC};;
		}
	}
	std::cout << "Метод Зейделя: превышен лимит итераций\n";
	return {};
}

		if (checkConvergence(x, x_old, eps)) {
			return x;
		}
	}
	throw std::runtime_error("Ошибка: превышен лимит итераций");
}
#endif //LAB_2__SEIDEL_METHOD_H_
