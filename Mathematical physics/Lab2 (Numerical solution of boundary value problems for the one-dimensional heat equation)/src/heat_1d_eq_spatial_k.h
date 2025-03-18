#ifndef HEAT_1D_EQ_SPATIAL_K_H
#define HEAT_1D_EQ_SPATIAL_K_H
#include <functional>
#include <variant>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h"

// template <typename T>
// using cond_type = std::variant<std::function<T(T)>, std::function<T(T, T)>, std::function<T(T, T, T)>>;

template <typename T>
std::vector<std::vector<T>> heat_eq_spatial(T rho, T c, T k1, T k2, T x1, T x2, T Q, T t0, T L, T u0,  std::function<T(T,T,T,T,T,T)> K, std::function<T(T,T,T)> initial_cond,  std::function<T(T)> left_cond, std::function<T(T)> right_cond, T tau, T h, T sigma)
{
	size_t num_points_in_space = std::round(L / h)+1;
	size_t num_points_over_time = std::round(t0 / tau)+1;
	std::vector<std::vector<T>> solution(1, std::vector<T>(num_points_in_space));
	std::vector<T> A(num_points_in_space-3);
	std::vector<T> a(num_points_in_space);
	std::vector<T> B(num_points_in_space-3);
	std::vector<T> C(num_points_in_space-2);
	std::vector<T> F(num_points_in_space-2);
	solution[0][0]=initial_cond(0 , u0, L);;
	for (size_t i = 1; i < num_points_in_space; i++){
		solution[0][i] = initial_cond(i*h , u0, L);
		a[i] = 2.0 * K((i - 1)*h, x1, x2, k1, k2, L) * K(i*h, x1, x2, k1, k2, L) / (K((i - 1)*h, x1, x2, k1, k2, L) + K(i*h, x1, x2, k1, k2, L));
		// a[i] = std::sqrt(K((i - 1)*h, x1, x2, k1, k2, L) * K(i*h, x1, x2, k1, k2, L));
	}
	for (size_t j = 1; j < num_points_over_time; j++)
	{
		for (size_t i = 1; i <= num_points_in_space-3; i++)
		{
			A[i-1] = (sigma/ h) * a[i];
			B[i-1] = (sigma/ h) * a[i+1];
		}
		for (size_t i = 1; i <= num_points_in_space-2; i++)
		{
			C[i-1] = + sigma/ h * a[i] + sigma/ h * a[i+1] + c*rho*h/tau;
			F[i-1] = c*rho*h/tau*solution[j-1][i]+(1-sigma)*(a[i + 1] * (solution[j-1][i + 1] - solution[j-1][i]) / h - a[i] * (solution[j-1][i] - solution[j-1][i-1]) / h);
		}
		// std::cout << A << std::endl << B << std::endl << C << std::endl << F << std::endl;
		// solution.emplace_back(TridiagonalMatrixAlgorithm(A,C,B,F));
		// solution[j][0] = left_cond(j*tau);
		std::vector<double> temp{left_cond(j * tau)};
		std::vector<double> internal_points = TridiagonalMatrixAlgorithm(A, C, B, F);
		temp.insert(temp.end(), internal_points.begin(), internal_points.end());
		// double P_tj1 = right_cond(tau*j, Q, t0);
		// double P_tj = right_cond(tau*(j-1), Q, t0);
		// double chi = (sigma * a[num_points_in_space-1] / h) / (c * rho * h / (2 * tau) + sigma * a[num_points_in_space-1] / h);
		// double numerator = c * rho * solution[j-1][num_points_in_space-1] * h / (2 * tau) + sigma * P_tj1 + (1 - sigma) * (P_tj - (solution[j-1][num_points_in_space-1] - solution[j-1][num_points_in_space-2]) / h);
		// double denominator = c * rho * h / (2 * tau) + sigma * a[num_points_in_space-1] / h;
		// double mu = numerator / denominator;
		// temp.push_back(chi*temp[num_points_in_space-2]+mu);
		temp.push_back(right_cond(j * tau));
		solution.emplace_back(std::move(temp));
	}
	return solution;
}
#endif //HEAT_1D_EQ_SPATIAL_K_H
