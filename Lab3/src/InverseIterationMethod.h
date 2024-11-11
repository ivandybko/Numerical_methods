#ifndef LAB3__INVERSEITERATIONMETHOD_H_
#define LAB3__INVERSEITERATIONMETHOD_H_
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include "../../Lab1/src/Gauss.h"
#include "../../Lab1/src/matrix_operations.h"
#include "../../Lab2/src/StoppingCriteria.h"

template <typename T>
std::vector<T> generateRandomVector(int size, int min, int max) {
	std::vector<T> vec(size);
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(min, max);

	for (auto& num : vec) {
		num = dist(gen);
	}
	return vec;
}

template <typename T>
bool isSymmetric(const std::vector<std::vector<T>>& matrix) {
	size_t n = matrix.size();
	for (const auto& row : matrix) {
		if (row.size() != n) {
			return false;
		}
	}
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			if (matrix[i][j] != matrix[j][i]) {
				return false;
			}
		}
	}
	return true;
}

template <typename T>
std::vector<std::unique_ptr<std::vector<T>>> InverseIterationMethod(const std::vector<std::vector<T>> &matrix, const std::vector<T> &eigenValues, T eps){
	size_t n{matrix.size()};
	std::vector<std::unique_ptr<std::vector<T>>> eigenVectors;
	for (auto eigenValue : eigenValues)
	{
		auto x{generateRandomVector<T>(n, 0, eigenValue+1)};
		std::vector<T> x_new(n,0);
		std::vector<std::vector<T>> A = matrix;
		for (size_t i = 0; i < n; ++i)
		{
			A[i][i]-=eigenValue;
		}
		while(true){
			x_new = Gauss(std::pair<std::vector<std::vector<T>>, std::vector<T>>(A,x));
			normalize(x_new);
			if (checkConvergence(x,x_new,eps) or checkConvergence(multiplyVectorByConstant<T>(x,-1),x_new,eps)){break;}
			x=x_new;
		}
		eigenVectors.push_back(std::make_unique<std::vector<T>>(x_new));
	}
	return eigenVectors;
}

template <typename T>
std::pair<std::vector<T>,std::vector<std::unique_ptr<std::vector<T>>>> InverseIterationMethod(const std::vector<std::vector<T>> &matrix, T eps){
	if(!isSymmetric(matrix)){
		throw std::invalid_argument("Матрица не симметрична. Отношение Рэлея не может быть использовано.");
	}
	size_t n{matrix.size()};
	std::vector<std::unique_ptr<std::vector<T>>> eigenVectors;
	std::vector<T> eigenValues;
	for (size_t i = 0; i < n; ++i)
	{
		auto x{generateRandomVector<T>(n, 0, 100)};
		std::vector<T> x_new(n,0);
		std::vector<std::vector<T>> A = matrix;
		normalize(x);
		T eigenValue = (A*x)*x;
		while(true){
			for (size_t j = 0; j < n; ++j)
			{
				A[j][j]-=eigenValue;
			}
			x_new = Gauss(std::pair<std::vector<std::vector<T>>, std::vector<T>>(A,x));
			for (size_t j = 0; j < n; ++j)
			{
				A[j][j]=matrix[j][j];
			}
			normalize(x_new);
			eigenValue = (A*x_new)*x_new;
			if (checkConvergence(x,x_new,eps) or checkConvergence(multiplyVectorByConstant<T>(x,-1),x_new,eps)){break;}
			x=x_new;
		}
		eigenValues.push_back(eigenValue);
		eigenVectors.push_back(std::make_unique<std::vector<T>>(std::move(x_new)));
	}
	return std::pair<std::vector<T>,std::vector<std::unique_ptr<std::vector<T>>>>(eigenValues,std::move(eigenVectors));
}

#endif //LAB3__INVERSEITERATIONMETHOD_H_
