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
std::vector<std::unique_ptr<std::vector<T>>> InverseIterationMethod(const std::vector<std::vector<T>> &matrix, const std::vector<T> &eigenValues, T eps){
	size_t n{eigenValues.size()};
	std::vector<std::unique_ptr<std::vector<T>>> eigenVectors;
	for (auto eigenValue : eigenValues)
	{
		auto x{generateRandomVector<T>(n, 0, eigenValue+1)};
		std::vector<T> x_new(n,0);
		std::vector<std::vector<T>> A = matrix;
		for (int i = 0; i < n; ++i)
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
#endif //LAB3__INVERSEITERATIONMETHOD_H_
