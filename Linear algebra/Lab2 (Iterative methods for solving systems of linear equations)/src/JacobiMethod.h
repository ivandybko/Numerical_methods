#ifndef LAB_2__JACOBIMETHOD_H_
#define LAB_2__JACOBIMETHOD_H_
#include <iostream>
#include <vector>
#include <utility>
#include "../../Lab1/src/matrix_operations.h"
#include "Solution.h"
#include "StoppingCriteria.h"

template <typename T>
Solution<T> JacobiMethod(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data, int max_iter =1000, T eps=1){
	auto *A = &data.first;
	auto *b = &data.second;
	size_t n = (*A).size();
	//Расчет С
	std::vector<std::vector<T>> C(n,std::vector<T>(n));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			C[i][j] = i!=j ? -(*A)[i][j]/(*A)[i][i] : 0;
		}
	}
	T norm_C{octahedralNorm(C)};
	//Нахождение решения
	std::vector<T> x(n, 0.0);
	for (size_t iter = 0; iter < max_iter; ++iter) {
		std::vector<T> x_old{x};
		T sum{0};
		for (int i = 0; i < n; ++i) {
			sum=0;
			for (int j = 0; j < n; ++j)
			{
				if(i != j){
					sum+=(*A)[i][j]*x_old[j];
				}
			}
			x[i]=((*b)[i]-sum)/(*A)[i][i];
		}
		if (checkConvergence(x, x_old, eps)){
			return Solution<T>{x, iter, norm_C};
		}
	}
	std::cout << "Метод Якоби: превышен лимит итераций\n";
	return {};
}
#endif //LAB_2__JACOBIMETHOD_H_
