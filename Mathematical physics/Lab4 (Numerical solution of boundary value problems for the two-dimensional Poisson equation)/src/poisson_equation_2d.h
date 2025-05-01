#ifndef POISSON_EQUATION_2D_H
#define POISSON_EQUATION_2D_H
#include <functional>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h"


template <typename T>
void printTridiagonalMatrix(const std::vector<T>& a,
						   const std::vector<T>& b,
						   const std::vector<T>& c) {
	size_t n = b.size(); // Размер матрицы n x n

	for (size_t i = 0; i < n; ++i) {
		std::cout << "{ ";
		for (size_t j = 0; j < n; ++j) {
			if (j == i) {
				// Главная диагональ
				std::cout << b[i] << ", ";
			} else if (j == i - 1 && i > 0) {
				// Поддиагональ (нижняя)
				std::cout << a[i - 1] << ", ";
			} else if (j == i + 1 && i < n - 1) {
				// Наддиагональ (верхняя)
				std::cout << c[i] << ", ";
			} else {
				// Вне диагоналей
				std::cout << "0, ";
			}
		}
		std::cout << " }," << std::endl;
	}
}

template <typename T>
struct BoundaryCondition {
	std::function<T(T)> func;
	bool is_flux;
};
 // Г1(x2=0) - bottom_bc, Г2(x2=Y) - top_bc, Г3(x1=0) - left_bc, Г4(x1=X) - right_bc
template <typename T>
std::vector<std::vector<std::vector<T>>> poisson_equation_2d(std::pair<T, T> upper_right_corner, T t0, std::function<T(T,T)> power_density, BoundaryCondition<T> left, BoundaryCondition<T> right, BoundaryCondition<T> bottom, BoundaryCondition<T> top, T tau, T h1 , T h2) {
	size_t num_points_in_space_1 = std::round(upper_right_corner.first / h1)+1;
	size_t num_points_in_space_2 = std::round(upper_right_corner.second / h2)+1;
	size_t num_points_over_time = std::round(t0 / tau)+1;
  	std::vector<std::vector<std::vector<T>>> solution(num_points_over_time, std::vector<std::vector<T>>(num_points_in_space_1, std::vector<T>(num_points_in_space_2, 0)));
	std::vector<T> A_1(num_points_in_space_1-1, 1.0/(h1*h1));
	std::vector<T> B_1(num_points_in_space_1-1, 1.0/(h1*h1));
	std::vector<T> C_1(num_points_in_space_1,   - 2.0/tau - 2.0/(h1*h1));
	std::vector<T> F_1(num_points_in_space_1);
	std::vector<T> A_2(num_points_in_space_2-1, 1.0/(h2*h2));
	std::vector<T> B_2(num_points_in_space_2-1, 1.0/(h2*h2));
	std::vector<T> C_2(num_points_in_space_2,   - 2.0/tau - 2.0/(h2*h2));
	std::vector<T> F_2(num_points_in_space_2);
	for (size_t j = 0; j < num_points_in_space_2; ++j)
	{
		if (left.is_flux){
			// solution[0][0][j] = (solution[t-1][0][j]/tau - 2*left.func(j*h2)/h1);
		}
		else{
			solution[0][0][j] = left.func(j*h2);
		}
		if (right.is_flux){

		}
		else{
			solution[0][num_points_in_space_1-1][j] = right.func(j*h2);
		}
	}
	for (size_t i = 0; i < num_points_in_space_1; ++i)
	{
		if (top.is_flux)
		{
			solution[0][i][num_points_in_space_2 - 1] = solution[0][i][num_points_in_space_2 - 2] + h2 * top.func(i*h1);
		}
		else
		{
			solution[0][i][num_points_in_space_2-1] = top.func(i*h1);
		}
		if (bottom.is_flux)
		{
			solution[0][i][0] = solution[0][i][1] - h2 * bottom.func(i*h1);

			// C_2[num_points_in_space_2-1] = - 2.0/(h2*h2);
			// A_2[num_points_in_space_2-2] = 2.0/(h2*h2);
			// F_2[num_points_in_space_2-1] = solution[t-1][i][num_points_in_space_2-1]/tau - 2*bottom.func(i*h1)/h2;
		}
		else
		{
			solution[0][i][0] = bottom.func(i*h1);
		}
	}
	for (size_t t = 1; t < num_points_over_time; t++)
	{
		std::vector<std::vector<T>> u_half(1, std::vector<T>(num_points_in_space_1, 0));
		for (size_t i = 0; i < num_points_in_space_1; i++)
		{
			if (bottom.is_flux)
			{
				// u_half[0][i] = (solution[t-1][i][0]*h2/tau + bottom.func(i*h1)) / (h2 / tau + 1.0/h2);
				u_half[0][i] = solution[t-1][i][0]/tau + bottom.func(i*h1)/h2;
			}
			else
			{
				u_half[0][i] = bottom.func(i*h1);
			}
		}
		for (size_t x2 = 1; x2 < num_points_in_space_2-1; x2++)
		{
			auto j=x2;
			if (left.is_flux)
			{
				// C_1[0] = - 2.0/(h1*h1);
				// C_1[0] = 2.0/(h1*h1);
				B_1[0] =  2.0/(h1*h1);
				// C_1[0] = -(1.0 / tau + 2.0 / (h1 * h1));
				F_1[0] = -(solution[t-1][0][j]/tau + 2*left.func(j*h2)/h1);
				// F_1[0] = 2*left.func(j*h2)/h1;
				// B_1[0] =0;
				// B_1[0] = 1.0/(h1*h1 /tau + 1);
				// F_1[0] = (solution[t-1][0][j] * h1 /tau + left.func(j*h2))/ (h1 / tau + 1.0/h1);
			}
			else
			{
				B_1[0] = 0;
				C_1[0] = 1;
				F_1[0] = left.func(j*h2);
				solution[t][0][j] = left.func(j*h2);
			}
			if (right.is_flux)
			{
				// C_1[num_points_in_space_1-1] = - 2.0/(h1*h1);
				A_1[num_points_in_space_1-2] = 2.0/(h1*h1);
				// C_1[num_points_in_space_1 - 1] = -(1.0 / tau + 2.0 / (h1 * h1));
				F_1[num_points_in_space_1-1] = -(solution[t-1][num_points_in_space_1-1][j]/tau + 2*right.func(j*h2)/h1);
				// A_1[num_points_in_space_1-2] = 0;
				// A_1[num_points_in_space_1-2] = 1.0/(h1*h1 /tau + 1);
				// F_1[num_points_in_space_1-1] = -(solution[t-1][num_points_in_space_1-1][j] * h1 /tau + right.func(j*h2))/ (h1 / tau + 1.0/h1);
			}
			else
			{
				A_1[num_points_in_space_1-2]=0;
				C_1[num_points_in_space_1-1]=1;
				F_1[num_points_in_space_1-1] = right.func(j*h2);
				solution[t][num_points_in_space_1-1][j] = right.func(j*h2);
			}
			for (size_t i = 1; i < F_1.size()-1; i++)
			{
				F_1[i] = - 2.0/tau * solution[t-1][i][x2] - (solution[t-1][i][x2-1] - 2*solution[t-1][i][x2] + solution[t-1][i][x2+1]) / (h2 * h2) - power_density(i * h1, x2 * h2);
			}
			// printTridiagonalMatrix(A_1, C_1,B_1);
			// std::cout << F_1 << std::endl;
			u_half.push_back(TridiagonalMatrixAlgorithm(A_1, C_1, B_1, F_1));
		}
		u_half.push_back(std::vector<T>(num_points_in_space_1, 0));
		for (size_t i = 0; i < num_points_in_space_1; i++)
		{
			if (top.is_flux)
			{
				// u_half[num_points_in_space_2-1][i] = (solution[t-1][i][num_points_in_space_2-1]*h2/tau + top.func(i*h1)) / (h2 / tau + 1.0/h2);
				u_half[num_points_in_space_2-1][i] = solution[t-1][i][num_points_in_space_2-1]/tau + 2*top.func(i*h1)/h2;
			}
			else
			{
				u_half[num_points_in_space_2-1][i] = top.func(i*h1);
			}
		}
		// std::cout << t << std::endl;
		// for (auto element : u_half)
		// {
		// 	std::cout << element << std::endl;
		// }
		transpose(u_half);
		for (size_t x1 = 1; x1 < num_points_in_space_1-1; x1++)
		{
			auto i=x1;
			if (bottom.is_flux) {
				B_2[0] = 2.0/(h2*h2);
				// C_2[0] = 1.0 / tau + 2.0 / (h2 * h2);
				// C_2[0] = -2.0/tau - 2.0/(h2*h2);
				F_2[0] = -(solution[t-1][i][0]/tau + 2*bottom.func(i*h1)/h2);
				// B_2[0] = 0;
				// B_2[0] =  -1.0/(h2*h2 / tau + 1);
				// F_2[0] = (solution[t-1][i][0]*h2/tau + bottom.func(i*h1)) / (h2 / tau + 1.0/h2);
			} else {
				B_2[0] = 0;
				C_2[0] = 1;
				F_2[0] = bottom.func(i*h1);
			}
			// Для x2 = Y (top)
			if (top.is_flux) {
				A_2[num_points_in_space_2-2] = 2.0/(h2*h2);
				// C_2[num_points_in_space_2-1] = 1.0/tau + 2.0/(h2*h2);
				// A_2[num_points_in_space_2-2] = 0;
				// A_2[num_points_in_space_2-2] = -1.0/(h2*h2 / tau + 1);
				// C_2[num_points_in_space_2-1] = -2.0/tau - 2.0/(h2*h2);
				F_2[num_points_in_space_2-1] = -(solution[t-1][i][num_points_in_space_2-1]/tau + 2*top.func(i*h1)/ h2);

				// F_2[num_points_in_space_2-1] = (solution[t-1][i][num_points_in_space_2-1]*h2/tau + top.func(i*h1)) / (h2 / tau + 1.0/h2);
			} else {
				A_2[num_points_in_space_2-2] = 0;
				C_2[num_points_in_space_2-1] = 1;
				F_2[num_points_in_space_2-1] = top.func(i*h1);
			}
			for (size_t j = 1; j < F_2.size()-1; j++)
			{
				F_2[j] = - 2.0/tau * u_half[x1][j] - (u_half[x1-1][j] - 2*u_half[x1][j] + u_half[x1+1][j]) / (h1 * h1) - power_density(x1 * h1, j * h2);
			}
			auto u = TridiagonalMatrixAlgorithm(A_2, C_2, B_2, F_2);
			for (size_t x2 = 0; x2 < num_points_in_space_2; x2++){
				solution[t][x1][x2] = u[x2];
			}
		}
		// if (left.is_flux){
		// 	// solution[t][0][0] = left.func(0);
		// 	// solution[t][0][num_points_in_space_2-1] = left.func(num_points_in_space_2*h2);
		// 	// for (size_t x2 = 0; x2 < num_points_in_space_2; x2++){
		// 	// 	solution[t][0][x2] = (solution[t-1][0][x2]/tau + 2*left.func(x2*h2)/h1);
		// 	// }
		// 	for (size_t x2 = 0; x2 < num_points_in_space_2; x2++){
		// 		// solution[t][num_points_in_space_1-1][x2] = (solution[t-1][num_points_in_space_1-1][x2]/tau - 2*right.func(x2*h2)/h1);
		// 		solution[t][0][x2] = u_half[0][x2];
		// 	}
		// }
		// else{
		// 	// solution[t][0][0] = left.func(0);
		// 	// solution[t][0][num_points_in_space_2-1] = left.func(num_points_in_space_2*h2);
		// 	for (size_t x2 = 0; x2 < num_points_in_space_2; x2++){
		// 		solution[t][0][x2] = left.func(x2*h2);
		// 	}
		// }
		// if (right.is_flux){
		// 	// solution[t][num_points_in_space_1-1][0] = right.func(0);
		// 	// solution[t][num_points_in_space_1-1][num_points_in_space_2-1] = right.func(num_points_in_space_2*h2);
		// 	for (size_t x2 = 0; x2 < num_points_in_space_2; x2++){
		// 		// solution[t][num_points_in_space_1-1][x2] = (solution[t-1][num_points_in_space_1-1][x2]/tau - 2*right.func(x2*h2)/h1);
		// 		solution[t][num_points_in_space_1-1][x2] = u_half[num_points_in_space_1-1][x2];
		// 	}
		// }
		// else{
		// 	solution[t][num_points_in_space_1-1][0] = right.func(0);
		// 	solution[t][num_points_in_space_1-1][num_points_in_space_2-1] = right.func(num_points_in_space_2*h2);
		// }
		for (size_t x2 = 0; x2 < num_points_in_space_2; x2++){
			// solution[t][num_points_in_space_1-1][x2] = (solution[t-1][num_points_in_space_1-1][x2]/tau - 2*right.func(x2*h2)/h1);
			solution[t][0][x2] = u_half[0][x2];
			solution[t][num_points_in_space_1-1][x2] = u_half[num_points_in_space_1-1][x2];
		}
		solution[t][0][0] = (solution[t][1][0] + solution[t][0][1]) - solution[t][1][1];
		solution[t][num_points_in_space_1-1][0] = (solution[t][num_points_in_space_1 - 2][0] + solution[t][num_points_in_space_1-1][1]) - solution[t][num_points_in_space_1-2][1];
		solution[t][0][num_points_in_space_2-1] = (solution[t][1][num_points_in_space_2-1] + solution[t][0][num_points_in_space_2-2]) - solution[t][1][num_points_in_space_2-2];
		solution[t][num_points_in_space_1-1][num_points_in_space_2-1] = solution[t][num_points_in_space_1-1][num_points_in_space_2-1] = (solution[t][num_points_in_space_1 - 2][num_points_in_space_2-1] + solution[t][num_points_in_space_1-1][num_points_in_space_2-1 - 1]) - solution[t][num_points_in_space_1-2][num_points_in_space_2-2];
		// transpose(solution[t]);

		// std::cout << "sol" << std::endl;
		// for (auto element : solution[t])
		// {
		// 	std::cout << element << std::endl;
		// }
	}
	return solution;
};
#endif //POISSON_EQUATION_2D_H
