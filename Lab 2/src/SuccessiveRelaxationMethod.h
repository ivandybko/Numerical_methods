#ifndef LAB_2_SRC_SUCCESSIVERELAXATIONMETHOD_H_
#define LAB_2_SRC_SUCCESSIVERELAXATIONMETHOD_H_
#include <iostream>
#include <vector>
#include <utility>
#include "Solution.h"
#include "StoppingCriteria.h"

template <typename T>
T calcNormC(const std::vector<std::vector<T>>& A, T omega) {
	size_t n = A.size();
	std::vector<std::vector<T>> L(n, std::vector<T>(n, 0));
	std::vector<std::vector<T>> D(n, std::vector<T>(n, 0));
	std::vector<std::vector<T>> U(n, std::vector<T>(n, 0));

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (i > j) {
				L[i][j] = A[i][j];
			} else if (i == j) {
				D[i][j] = A[i][j];
			} else {
				U[i][j] = A[i][j];
			}
		}
	}

	std::vector<std::vector<T>> InvDL(n, std::vector<T>(n, 0));
	for (size_t i = 0; i < n; ++i) {
		InvDL[i][i] = 1 / (D[i][i] + omega * L[i][i]);
		for (size_t j = 0; j < i; ++j) {
			T sum = 0;
			for (size_t k = j; k < i; ++k) {
				sum += (omega * L[i][k]) * InvDL[k][j];
			}
			InvDL[i][j] = -sum / (D[i][i] + omega * L[i][i]);
		}
	}

	std::vector<std::vector<T>> DU(n, std::vector<T>(n, 0));
	for (size_t i = 0; i < n; ++i) {
		DU[i][i] = (1 - omega) * D[i][i];
		for (size_t j = i + 1; j < n; ++j) {
			DU[i][j] = -omega * U[i][j];
		}
	}

	std::vector<std::vector<T>> C(n, std::vector<T>(n, 0));
	multiply(InvDL, DU, C);
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
Solution<T> SuccessiveRelaxationMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter =1000, T eps=1e-4, T omega=1.1){
	auto *A = &data.first;
	auto *b = &data.second;
	size_t n = (*A).size();
	std::vector<T> x(n, 0.0);
	T normC{calcNormC<T>(*A,omega)};
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
			x[i] =(1-omega)*x_old[i]+omega*((*b)[i] - sum) / (*A)[i][i];
		}
		if (checkConvergence(x, x_old, eps)) {
			return Solution<T>{x, iter , normC};;
		}
	}
	std::cout << "Метод релаксации: превышен лимит итераций\n";
	return {};
}
#endif //LAB_2_SRC_SUCCESSIVERELAXATIONMETHOD_H_
