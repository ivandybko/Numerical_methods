#ifndef LAB3__HESSENBERGREDUCTION_H_
#define LAB3__HESSENBERGREDUCTION_H_
#include <vector>
#include <iostream>
#include <cmath>

template <typename T>
void givensRotation(T a, T b, T& c, T& s) {
	T r = std::sqrt(a * a + b * b);
	c = a / r;
	s = -b / r;
}

template <typename T>
void applyGivensLeft(std::vector<std::vector<T>>& A, T c, T s, size_t i, size_t j) {
	size_t n = A.size();
	for (size_t k = 0; k < n; ++k) {
		T temp = c * A[i][k] - s * A[j][k];
		A[j][k] = s * A[i][k] + c * A[j][k];
		A[i][k] = temp;
	}
}

template <typename T>
void applyGivensRight(std::vector<std::vector<T>>& A, T c, T s, size_t i, size_t j) {
	size_t n = A.size();
	for (size_t k = 0; k < n; ++k) {
		T temp = c * A[k][i] - s * A[k][j];
		A[k][j] = s * A[k][i] + c * A[k][j];
		A[k][i] = temp;
	}
}

template <typename T>
void reduceToHessenbergWithGivens(std::vector<std::vector<T>>& A) {
	size_t n = A.size();
	for (size_t j = 0; j < n - 2; ++j) {
		for (size_t i = n - 1; i > j + 1; --i) {
//			if (A[i][j] != 0) {
				T c, s;
				givensRotation(A[i - 1][j], A[i][j], c, s);
				applyGivensLeft(A, c, s, i - 1, i);
				applyGivensRight(A, c, s, i - 1, i);
//			}
		}
	}
}



#endif //LAB3__HESSENBERGREDUCTION_H_
