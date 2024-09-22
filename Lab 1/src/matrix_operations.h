#ifndef LAB_1_SRC_MATRIX_OPERATIONS_H_
#define LAB_1_SRC_MATRIX_OPERATIONS_H_
#include <iostream>
#include <vector>
#include <utility>
template <typename T>
std::vector<std::vector<T>>* inverse(const std::vector<std::vector<T>>& A){
	size_t n = A.size();
	std::vector<std::vector<T>> identity(n, std::vector<T>(n, 0));
	for (size_t i = 0; i < n; ++i) {
		identity[i][i] = 1;
	}
	std::vector<T> solution(n);
	auto* inverse = new std::vector<std::vector<T>>(n, std::vector<T>(n));
	for (size_t i = 0; i < n; ++i) {
		std::pair<std::vector<std::vector<T>>, std::vector<T>> data(A, identity[i]);
		solution = Gauss(data);
		for (size_t j = 0; j < n; ++j) {
			(*inverse)[j][i] = solution[j];
		}
	}
	return inverse;
}

template <typename T>
void transpose(std::vector<std::vector<T>>& A) {
	size_t n = A.size();
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			std::swap(A[i][j], A[j][i]);
		}
	}
}
template <typename T>
void multiply(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B, std::vector<std::vector<T>>& res) {
	size_t n = A.size();
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			for (int k = 0; k < n; ++k) {
				res[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

template <typename T>
T octahedralNorm(const std::vector<std::vector<T>>& A) {
	size_t n = A.size();
	if (n == 0) return 0;
	T maxColSum = 0;
	for (size_t j = 0; j < n; ++j) {
		T colSum = 0;
		for (size_t i = 0; i < n; ++i) {
			colSum += std::abs(A[i][j]);
		}
		if (colSum > maxColSum) {
			maxColSum = colSum;
		}
	}
	return maxColSum;
}

template <typename T>
T cubicNorm(const std::vector<std::vector<T>>& A) {
	size_t n = A.size();
	if (n == 0) return 0;
	T maxRowSum = 0;
	for (size_t i = 0; i < n; ++i) {
		T rowSum = 0;
		for (size_t j = 0; j < n; ++j) {
			rowSum += std::abs(A[i][j]);
		}
		if (rowSum > maxRowSum) {
			maxRowSum = rowSum;
		}
	}
	return maxRowSum;
}
#endif //LAB_1_SRC_MATRIX_OPERATIONS_H_
