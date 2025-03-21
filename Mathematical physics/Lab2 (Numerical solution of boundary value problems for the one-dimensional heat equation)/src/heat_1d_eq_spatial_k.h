#ifndef HEAT_1D_EQ_SPATIAL_K_H
#define HEAT_1D_EQ_SPATIAL_K_H
#include <functional>
#include <variant>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h"

// template <typename T>
// using cond_type = std::variant<std::function<T(T)>, std::function<T(T, T)>, std::function<T(T, T, T)>>;

template <typename T>
std::vector<std::vector<T>> heat_eq_spatial(T rho, T c, T k1, T k2, T x1, T x2, T Q, T t0, T L, T u0,  std::function<T(std::vector<T>)> K, std::function<T(std::vector<T>)> initial_cond,  std::function<T(std::vector<T>)> left_cond, bool stream_left, std::function<T(std::vector<T>)> right_cond, bool stream_right, T tau, T h, T sigma)
{
	size_t num_points_in_space = std::round(L / h)+1;
	size_t num_points_over_time = std::round(t0 / tau)+1;
	std::vector<std::vector<T>> solution(1, std::vector<T>(num_points_in_space));
	std::vector<T> A(num_points_in_space-3);
	std::vector<T> a(num_points_in_space);
	std::vector<T> B(num_points_in_space-3);
	std::vector<T> C(num_points_in_space-2);
	std::vector<T> F(num_points_in_space-2);
	solution[0][0]=initial_cond({0 , L, u0});;
	for (size_t i = 1; i < num_points_in_space; i++){
		solution[0][i] = initial_cond({i*h , L, u0});
		// a[i] = 2.0 * K({(i - 1)*h, L, x1, x2, k1, k2, L}) * K({i*h,L, x1, x2, k1, k2}) / (K({(i - 1)*h, L, x1, x2, k1, k2, L}) + K({i*h, L, x1, x2, k1, k2}));
		a[i] = std::sqrt(K({(i - 1)*h, L, x1, x2, k1, k2, L}) * K({i*h, L, x1, x2, k1, k2, L}));
		// a[i] = std::sqrt(K((i - 1)*h, x1, x2, k1, k2, L) * K(i*h, x1, x2, k1, k2, L));
	}
	for (size_t j = 1; j < num_points_over_time; j++)
	{
		for (size_t i = 1; i <= num_points_in_space-3; i++)
		{
			A[i-1] = (sigma/ h) * a[i+1];
			B[i-1] = (sigma/ h) * a[i+1];
		}
		for (size_t i = 1; i <= num_points_in_space-2; i++)
		{
			C[i-1] = - (sigma/ h * a[i] + sigma/ h * a[i+1] +  c*rho*h/tau);
			F[i-1] = - (c*rho*h/tau*solution[j-1][i] + (1-sigma)*(a[i + 1] * (solution[j-1][i + 1] - solution[j-1][i]) / h - a[i] * (solution[j-1][i] - solution[j-1][i-1]) / h));
		}
		std::vector<T> internal_points = TridiagonalMatrixAlgorithm(A, C, B, F);
		std::vector<T> temp{};
		// if (stream_left)
		// {
		// 	double P_tj1 = left_cond({tau*j, Q, t0});
		// 	double P_tj = left_cond({tau*(j-1), Q, t0});
		// 	double chi = (sigma * a[1] / h) / (c * rho * h / (2 * tau) + sigma * a[1] / h);
		// 	double numerator = c * rho * solution[j-1][0] * h / (2 * tau) + sigma * P_tj1 + (1 - sigma) * (P_tj - (solution[j-1][1] - solution[j-1][0]) / h);
		// 	double denominator = c * rho * h / (2 * tau) + sigma * a[1] / h;
		// 	double mu = numerator / denominator;
		// 	temp.push_back(chi*internal_points[0]+mu);
		// }
		if (stream_left)
		{
			T P_tj1 = left_cond({tau*j, Q, t0});
			T P_tj = left_cond({tau*(j-1), Q, t0});
			T left= ((h*h+ 2*(-1+sigma)*tau*a[1])*solution[j-1][0]+2*tau*(h*(P_tj-P_tj*sigma+P_tj1* sigma)-(-1+sigma)*a[1]*solution[j-1][1]*sigma*a[1]*internal_points[0]))/(h*h+2*sigma*tau * a[1]);
			temp.push_back(left);
		}

		// solution[j][0] = left_cond(j*tau);
		else
		{
			temp.push_back(left_cond({j * tau}));
		}
		temp.insert(temp.end(), internal_points.begin(), internal_points.end());
		if (stream_right)
		{
			T P_tj1 = right_cond({tau*j, Q, t0});
			T P_tj = right_cond({tau*(j-1), Q, t0});
			T chi = (sigma * a[num_points_in_space-1] / h) / (c * rho * h / (2 * tau) + sigma * a[num_points_in_space-1] / h);
			T numerator = c * rho * solution[j-1][num_points_in_space-1] * h / (2 * tau) + sigma * P_tj1 + (1 - sigma) * (P_tj - (solution[j-1][num_points_in_space-1] - solution[j-1][num_points_in_space-2]) / h);
			T denominator = c * rho * h / (2 * tau) + sigma * a[num_points_in_space-1] / h;
			T mu = numerator / denominator;
			temp.push_back(chi*temp[num_points_in_space-2]+mu);
		}
		else{
			temp.push_back(right_cond({j * tau}));
		}
		solution.emplace_back(std::move(temp));
	}
	return solution;
}
#endif //HEAT_1D_EQ_SPATIAL_K_H
