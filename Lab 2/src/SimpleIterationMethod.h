#ifndef LAB_2__SIMPLEITERATIONMETHOD_H_
#define LAB_2__SIMPLEITERATIONMETHOD_H_
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "../../Lab 1/src/matrix_operations.h"
template <typename T>
bool checkConvergence(const std::vector<T> & b_minus_Ax, const std::vector<T> &x, const std::vector<T> &x_old, T eps)
{
	size_t n = x.size();
	T criteria1{ 0.0 }, criteria2{ 0.0 };
	for (size_t i = 0; i < n; ++i)
	{
		criteria1 += abs(b_minus_Ax[i]);
		criteria2 += abs(x[i] - x_old[i]);
	}
	return criteria1 < eps or criteria2 < eps;
}

template <typename T>
std::vector<T> SimpleIterationMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter = 1000, T tau = 0.001, T eps =0.0001){
	auto A = data.first;
	auto b = data.second;
	size_t n = A.size();
	auto norm_C{std::trunc(octahedralNorm(-tau*A+ identityMatrix<T>(n))* 1e3)/ 1e3};
	if (norm_C  > 1){
		std :: cout << "Норма C: " << norm_C << '\n';
		throw std::invalid_argument("Параметр τ выбран неверно. Сходимость невозможна.");
	}
	std::vector<T> x(n, 0.0);
	for (int iter = 0; iter < max_iter; ++iter) {
		std::vector<double> x_old{x};
		auto b_minus_Ax = b - A*x_old;
		x= tau*b_minus_Ax + x_old;
		T criteria1{0.0}, criteria2{0.0};
		if (checkConvergence<T>(b_minus_Ax, x, x_old, eps)){return x;}
	}
	throw std::runtime_error("Ошибка: превышен лимит итераций");
}
#endif //LAB_2__SIMPLEITERATIONMETHOD_H_
