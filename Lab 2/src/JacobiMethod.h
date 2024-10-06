#ifndef LAB_2__JACOBIMETHOD_H_
#define LAB_2__JACOBIMETHOD_H_
#include <iostream>
#include <vector>
#include <utility>
#include "../../Lab 1/src/matrix_operations.h"
//template <typename T>
//bool checkConvergence(const std::vector<T> & b_minus_Ax, const std::vector<T> &x, const std::vector<T> &x_old, T eps)
//{
//	size_t n = x.size();
//	T criteria1{ 0.0 }, criteria2{ 0.0 };
//	for (size_t i = 0; i < n; ++i)
//	{
//		criteria1 += abs(b_minus_Ax[i]);
//		criteria2 += abs(x[i] - x_old[i]);
//	}
//	return criteria1 < eps or criteria2 < eps;
//}

template <typename T>
std::vector<T> JacobiMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter =1000, T eps=1e-4){
	auto A = data.first;
	auto b = data.second;
	size_t n = A.size();
	std::vector<T> x(n, 0.0);
	std::vector<T> D(n, 0.0);
	for (int i = 0; i < n; ++i)
	{
		D[i]=A[i][i];
		A[i][i]=0;
	}
	for (int iter = 0; iter < max_iter; ++iter) {
		std::vector<T> x_old{x};
		auto b_minus_Ax = (b-(A*x_old)); // Здесь A=L+U
		std::vector<T> D_dot_x_old(n,0);
		for (int i = 0; i < n; ++i) {
			x[i] = b_minus_Ax[i] / D[i];
			D_dot_x_old[i] = D[i]*x_old[i];
		}
		if (checkConvergence(b_minus_Ax-D_dot_x_old,x, x_old, eps)) {
			return x;
		}
	}
	throw std::runtime_error("Ошибка: превышен лимит итераций");
}
#endif //LAB_2__JACOBIMETHOD_H_
