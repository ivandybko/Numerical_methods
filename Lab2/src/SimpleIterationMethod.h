#ifndef LAB_2__SIMPLEITERATIONMETHOD_H_
#define LAB_2__SIMPLEITERATIONMETHOD_H_
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include "../../Lab 1/src/matrix_operations.h"
#include "Solution.h"
#include "StoppingCriteria.h"

template <typename T>
Solution<T> SimpleIterationMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter = 1000, T tau = 0.001, T eps =0.0001){
	auto *A = &data.first;
	auto *b = &data.second;
	size_t n = (*A).size();
	auto norm_C{octahedralNorm(-tau*(*A) + identityMatrix<T>(n))};
//	if (norm_C  >= 1){
//		std :: cout << "Норма C: " << norm_C << '\n';
//		throw std::invalid_argument("Параметр τ выбран неверно. Сходимость невозможна.");
//	}
	std::vector<T> x(n, 0.0);
	for (size_t iter = 0; iter < max_iter; ++iter) {
		std::vector<double> x_old{x};
		auto b_minus_Ax = (*b) - (*A)*x_old;
		x= tau*b_minus_Ax + x_old;
		if (checkConvergence<T>(b_minus_Ax, eps) or checkConvergence<T>(x, x_old, eps,norm_C)){return Solution<T>{x, iter ,norm_C};}
	}
	std::cout << "Метод простой итерации: превышен лимит итераций\n";
	return {};
}
#endif //LAB_2__SIMPLEITERATIONMETHOD_H_