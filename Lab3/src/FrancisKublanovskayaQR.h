#ifndef LAB3__FRANCISKUBLANOVSKAYAQR_H_
#define LAB3__FRANCISKUBLANOVSKAYAQR_H_
#include <iostream>
#include <vector>
#include "../../Lab1/src/QR decomposition.h"
#include "../../Lab1/src/matrix_operations.h"
#include "HessenbergReduction.h"

template <typename T>
std::unique_ptr<std::vector<std::vector<T>>> multiplyRQ(const std::vector<std::vector<T>>& R, const std::vector<std::vector<T>>& Q) {
	size_t n = R.size();
	auto result = std::make_unique<std::vector<std::vector<T>>>(n, std::vector<T>(n, 0));
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i; j < n; ++j) {
			for (size_t k = 0; k < n; ++k) {
				(*result)[i][k] += R[i][j] * Q[j][k];
			}
		}
	}
	return result;
}

template <typename T>
std::vector<T> QRAlgorithm(const std::vector<std::vector<T>> &matrix, T eps, bool ShiftMode=false, bool HessenbergReduction=false){
	T shift{0};	int n = matrix.size(); std::vector<std::vector<T>> A{matrix};	std::vector<T> eigenvalues;
	if (HessenbergReduction){
		reduceToHessenbergWithGivens(A);
	}
	while (n > 0){
		if (ShiftMode){
			shift = A[n-1][n-1];
			for (int i = 0; i < n; ++i)
			{
				A[i][i]-=shift ;
			}
		}
		auto QR = QRdecomposition<T>(A);
		auto A_new = multiplyRQ((*QR).second, (*QR).first);
		T last_row{0};
		for (int i = 0; i < n; ++i)
		{
			(*A_new)[i][i] += shift;
			last_row+=(*A_new)[n-1][i];
		}
		if (abs(last_row - (*A_new)[n-1][n-1])  < eps){
			eigenvalues.push_back((*A_new)[n-1][n-1]);
			(*A_new).pop_back();
			for (auto& row : *A_new) {
				row.resize(row.size() - 1);
			}
			n-=1;
		}
		A = *A_new;
	}
   return eigenvalues;
}
#endif //LAB3__FRANCISKUBLANOVSKAYAQR_H_
