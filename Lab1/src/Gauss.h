#ifndef LAB_1_SRC_GAUSS_H_
#define LAB_1_SRC_GAUSS_H_
#include <iostream>
#include <vector>
#include <utility>
#include <limits> 
template <typename T>
std::vector<T> Gauss(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data){
	auto A = data.first;
	auto b = data.second;
	size_t n = A.size();
	std::vector<T> y(n);
	std::vector<T> x(n);
	std::vector<std::vector<T>> C(n, std::vector<T>(n));
	for (size_t i = 0; i < n; ++i)
	{
		if (A[i][i] == 0 ){
			T max{-std::numeric_limits<T>::max()};
			size_t max_k = i;
			for (size_t k=i+1; k < n; ++k){
				if (A[i][k] > max){max=A[i][k]; max_k=k;}
			}
			if (max_k == i){
				return std::vector<T>();
			}
			std::swap(A[i], A[max_k]);
			std::swap(b[i], b[max_k]);
		}
		y[i]=b[i]/A[i][i];
		for (size_t j = i+1; j < n ; ++j)
		{
			C[i][j]=A[i][j]/A[i][i];
		}
		for (size_t i1 = i+1; i1 < n; ++i1)
		{
			b[i1]-=A[i1][i]*y[i];
			for (size_t j1 = i+1; j1 < n; ++j1)
			{
				A[i1][j1]-=A[i1][i]*C[i][j1];
			}
		}
	}
	x[n-1]=b[n-1]/A[n-1][n-1];
//	for (size_t i = n-1; i-- > 0; )
//	{
////		if (C[i][i]==0){
////			return std::vector<T>();
////		}
//	}
	for (size_t i = n-1; i-- > 0; )
	{
		T sum{0};
		for (size_t j = i; j < n; ++j)
		{
			sum+=C[i][j]*x[j];
		}
		x[i]=y[i]-sum;
	}
	return x;
}
#endif //LAB_1_SRC_GAUSS_H_
