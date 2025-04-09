#ifndef HEAT_1D_EQ_SPATIAL_K_H
#define HEAT_1D_EQ_SPATIAL_K_H
#include <functional>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h"

template <typename T>
std::vector<std::vector<T>> heat_eq_spatial(T rho, T c, T k1, T k2, T x1, T x2, T Q, T t0, T L, T u0,  std::function<T(std::vector<T>)> K, std::function<T(std::vector<T>)> initial_cond,  std::function<T(std::vector<T>)> left_cond, bool stream_left, std::function<T(std::vector<T>)> right_cond, bool stream_right, T tau, T h, T sigma)
{
	size_t num_points_in_space = std::round(L / h)+1;
	size_t num_points_over_time = std::round(t0 / tau)+1;
	std::vector<std::vector<T>> solution(1, std::vector<T>(num_points_in_space));
	std::vector<T> A(num_points_in_space-1);
	std::vector<T> a(num_points_in_space+1);
	std::vector<T> B(num_points_in_space-1);
	std::vector<T> C(num_points_in_space);
	std::vector<T> F(num_points_in_space);
	solution[0][0]=initial_cond({0 , L, u0});;
	for (size_t i = 0; i < num_points_in_space; i++){
		solution[0][i] = initial_cond({i*h , L, u0});
		a[i] = 2.0 * K({(i - 1)*h, L, x1, x2, k1, k2, L}) * K({i*h,L, x1, x2, k1, k2}) / (K({(i - 1)*h, L, x1, x2, k1, k2, L}) + K({i*h, L, x1, x2, k1, k2}));
	}
	for (size_t j = 1; j < num_points_over_time; j++)
	{
		for (size_t i = 1; i < num_points_in_space-1; i++)
		{
			A[i-1] = (sigma/ h) * a[i];
			B[i] = (sigma/ h) * a[i+1];
		}
		for (size_t i = 1; i < num_points_in_space-1; i++)
		{
			C[i] = - ((sigma/ h) * a[i] + (sigma/ h) * a[i+1] +  c*rho*h/tau);
			F[i] = - (c*rho*h/tau*solution[j-1][i] + (1-sigma)*(a[i + 1] * (solution[j-1][i + 1] - solution[j-1][i]) / h - a[i] * (solution[j-1][i] - solution[j-1][i-1]) / h));
		}
		if (stream_left)
		{
			C[0] = c * rho * h / 2 + tau * sigma * K({h/2,L, x1, x2, k1, k2}) / h;
			B[0] = -tau * sigma * K({h/2,L, x1, x2, k1, k2}) / h;
			F[0] = c * rho * h / 2 * solution[j-1][0]  + tau * (1 - sigma) * (K({h/2,L, x1, x2, k1, k2}) * (solution[j-1][1]  - solution[j-1][0] ) / h - left_cond({tau*j, Q, t0})) + tau * sigma * left_cond({tau*(j+1), Q, t0});
		}
		else
		{
			B[0]=0;
			C[0]=1;
			F[0]=left_cond({tau*j, Q, t0});
		}
		if (stream_right)
		{
			A[num_points_in_space-2] = -(sigma * a[num_points_in_space-1] / h) / ((c * rho * h) / (2 * tau) + sigma * a[num_points_in_space-1] / h);
			C[num_points_in_space-1] = 1;
			F[num_points_in_space-1] = ((c * rho * solution[j-1][num_points_in_space-1] * h) / (2 * tau) + sigma * right_cond({tau*(j+1), Q, t0}) + (1 - sigma) * (right_cond({tau*j, Q, t0}) - a[num_points_in_space-1] * (solution[j-1][num_points_in_space-1] - solution[j-1][num_points_in_space-2]) / h)) /
				(c * rho * h / (2 * tau) + sigma * a[num_points_in_space-1] / h);
		}
		else{
			A[num_points_in_space-2]=0;
			C[num_points_in_space-1]=1;
			F[num_points_in_space-1]=right_cond({tau*j, Q, t0});
		}
		std::vector<T> internal_points = TridiagonalMatrixAlgorithm(A, C, B, F);
 		solution.emplace_back(std::move(internal_points));
	}
	return solution;
}

template <typename T>
std::vector<std::vector<T>> heat_eq_nonlinear(T rho, T c, T alpha, T beta, T gamma, T t0, T L, T u0,  std::function<T(std::vector<T>)> K, std::function<T(std::vector<T>)> initial_cond,  std::function<T(std::vector<T>)> left_cond, bool stream_left, std::function<T(std::vector<T>)> right_cond, bool stream_right, T tau, T h, T max_iter=1000, T epsilon=1e-6)
{
	size_t num_points_in_space = std::round(L / h);
	size_t num_points_over_time = std::round(t0 / tau);
	std::vector<std::vector<T>> solution(1, std::vector<T>(num_points_in_space+1));
	std::vector<T> A(num_points_in_space);
	std::vector<T> a(num_points_in_space+1);
	std::vector<T> B(num_points_in_space+1);
	std::vector<T> C(num_points_in_space);
	std::vector<T> F(num_points_in_space+1);
	solution[0][0]=initial_cond({0 , L, u0});;
	std::vector<T> u_prev(num_points_in_space + 1, 0.0), u_current(num_points_in_space + 1);
	for (size_t i = 1; i < num_points_in_space; i++){
		solution[0][i] = initial_cond({i*h , L, u0});
	}
	for (int n = 1; n < num_points_over_time; ++n)
	{
		T t = n * tau;
		for (int iter = 0; iter < max_iter; ++iter)
		{
			u_prev = u_current;
			for (int i = 1; i < num_points_in_space; ++i) {
				T Ki = K({u_prev[i], alpha, beta, gamma});
				T Ki1 = K({u_prev[i + 1], alpha, beta, gamma});
				T Ki_1 = K({u_prev[i - 1], alpha, beta, gamma});
				T a_i = 0.5 * (Ki + Ki_1);
				T a_i1 = 0.5 * (Ki1 + Ki);
				A[i-1] = a_i / h;
				C[i] = a_i1 / h;
				B[i] = -A[i-1] - C[i] - c * rho *h / tau;
				F[i] = -c * rho * h * u_prev[i] / tau;
				if (i * h >= c * t)
				{
					A[i-1] = 0;
					B[i] = 1;
					C[i] = 0;
					F[i] = right_cond({});
				}
			}
			B[0] = 1;
			C[0] = 0;
			F[0] = left_cond({t, beta, gamma, c});
			A[num_points_in_space-1] = 0;
			B[num_points_in_space] = 1;
			F[num_points_in_space] = right_cond({});
			u_current = TridiagonalMatrixAlgorithm(A, B, C, F);
			double max_diff = 0.0;
			for (int i = 0; i <= num_points_in_space; ++i) {
				max_diff = std::max(max_diff, abs(u_current[i] - u_prev[i]));
			}
			if (max_diff < epsilon) break;
		}
		for (int i = 0; i <= num_points_in_space; i++)
		{
			u_prev[i] = u_current[i];
		}
		solution.emplace_back(u_current);
	}
	return solution;
}
#endif //HEAT_1D_EQ_SPATIAL_K_H
