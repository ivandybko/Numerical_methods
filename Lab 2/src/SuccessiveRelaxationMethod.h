#ifndef LAB_2_SRC_SUCCESSIVERELAXATIONMETHOD_H_
#define LAB_2_SRC_SUCCESSIVERELAXATIONMETHOD_H_
#include <iostream>
#include <vector>
#include <utility>

template <typename T>
std::vector<T> SuccessiveRelaxationMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter =1000, T eps=1e-4, T omega=1.1){
	auto A = data.first;
	auto b = data.second;
	size_t n = A.size();
	std::vector<T> x(n, 0.0);
	for (int iter = 0; iter < max_iter; ++iter) {
		std::vector<T> x_old{x};
		for (size_t i = 0; i < n; ++i) {
			T sum{0.0};
			for (size_t j = 0; j < i; ++j) {
				sum += A[i][j] * x[j];
			}
			for (size_t j = i + 1; j < n; ++j) {
				sum += A[i][j] * x_old[j];
			}
			x[i] =(1-omega)*x_old[i]+omega*(b[i] - sum) / A[i][i];
		}
		if (checkConvergence(x, x_old, eps)) {
			std::cout << "Iter: " << iter;
			return x;
		}
	}
	throw std::runtime_error("Ошибка: превышен лимит итераций");
}
#endif //LAB_2_SRC_SUCCESSIVERELAXATIONMETHOD_H_
