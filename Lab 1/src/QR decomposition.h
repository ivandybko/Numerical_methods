#ifndef LAB_1_SRC_QR_DECOMPOSITION_H_
#define LAB_1_SRC_QR_DECOMPOSITION_H_
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
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
void rotate(std::vector<std::vector<T>>& transit, std::vector<std::vector<T>>& A, int i, int j) {
	T c = A[i][i]/std::sqrt(A[i][i]*A[i][i]+A[j][i]*A[j][i]);
	T s = A[j][i]/std::sqrt(A[i][i]*A[i][i]+A[j][i]*A[j][i]);
	size_t n = A.size();
	for (int col = 0; col < n; ++col) {
		T temp = c * A[i][col] + s * A[j][col];
		A[j][col] = -s * A[i][col] + c * A[j][col];
		A[i][col] = temp;
		T temp1=c * transit[i][col] + s * transit[j][col];
		T temp2=(-s) * transit[i][col] + c * transit[j][col];
		transit[i][col] = temp1;
		transit[j][col] = temp2;
	}
}

template <typename T>
std::vector<T> QRdecompostition(const std::pair<std::vector<std::vector<T>>, std::vector<T>> &data){
	auto A = data.first;
	auto b = data.second;
	size_t n = A.size();
	std::vector<std::vector<T>> transit(n,std::vector<T>(n,0));
	for (int i = 0; i < n; ++i)
	{
		transit[i][i]=1;
	}
	for (int i = 0; i < n-1; ++i)
	{
		for (int j = i+1; j < n; ++j)
		{
			if (A[j][i] != 0 or A[i][i]!=0){
				rotate(transit,A,i,j);
			}
		}
	}
	std::vector<T> x(n);
	std::vector<T> b1(n,0);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			b1[i]+=transit[i][j]*b[j];
		}
	}
	x[n-1]=b1[n-1]/A[n-1][n-1];
	for (size_t i = n-1; i-- > 0; )
	{
		T sum{0};
		for (size_t j = i; j < n; ++j)
		{
			sum+=A[i][j]*x[j];
		}
		x[i]=(b1[i]-sum)/A[i][i];
	}
	transpose(transit);
//	std::vector<std::vector<T>> test(n,std::vector<T>(n,0));
//	multiply(transit,A,test);
//	std::cout << "Q:\n";
//	for (auto row : transit)
//	{
//		for (auto element : row)
//		{
//			std::cout << element << ' ';
//		}
//		std::cout << '\n';
//	}
//	std::cout << "R:\n";
//	for (auto row : A)
//	{
//		for (auto element : row)
//		{
//			std::cout << element << ' ';
//		}
//		std::cout << '\n';
//	}
//	std::cout << "QR:\n";
//	for (auto row : test)
//	{
//		for (auto element : row)
//		{
//			std::cout << element << ' ';
//		}
//		std::cout << '\n';
//	}
	return x;
}
#endif //LAB_1_SRC_QR_DECOMPOSITION_H_
