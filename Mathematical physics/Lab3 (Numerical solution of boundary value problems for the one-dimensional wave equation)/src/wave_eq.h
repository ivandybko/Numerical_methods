#ifndef WAVE_EQ_H
#define WAVE_EQ_H

#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"

template <typename T>
std::vector<std::vector<T>> wave_eq(T a, T t0, T left, T right, std::function<T(T)> initial_cond, std::function<T(T)> initial_speed, std::function<T(T)> boundary_left, std::function<T(T)> boundary_right , T tau, T h, std::function<T(T)> initial_cond_d = nullptr)
{
	size_t num_points_in_space = std::round((right-left) / h)+1;
	size_t num_points_over_time = std::round(t0 / tau)+1;
	std::vector<std::vector<T>> solution(num_points_over_time, std::vector<T>(num_points_in_space));
	for (size_t i = 0; i < 2; i++)
	{
		solution[i][0] = boundary_left(i*tau);
		solution[i][num_points_in_space-1] = boundary_right(i*tau);
	}
	if (initial_cond_d)
	{
		T point{left+h};
		for (int j = 1; j < num_points_in_space-1; j++){
			solution[0][j] = initial_cond(point);
			solution[1][j] = solution[0][j] + tau * initial_speed(point) + ( a * a * tau * tau) / 2 * initial_cond_d(point);
			point+=h;
		}
	}
	else
	{
		T point{left+h};
		for (size_t j = 1; j < num_points_in_space-1; j++){
			solution[0][j] = initial_cond(point);
			solution[1][j] = solution[0][j] + tau * initial_speed(point) + ( a * a * tau * tau) / 2 * (initial_cond(point+h) -2*initial_cond(point) + initial_cond(point-h)) / ( h * h);
			point+=h;
		}
	}
	for (size_t i = 2; i < num_points_over_time; i++)
	{
		solution[i][0] = boundary_left(i*tau);
		solution[i][num_points_in_space-1] = boundary_right(i*tau);
		for (size_t j = 1; j < num_points_in_space-1; j++)
		{
			solution[i][j] = (a * a * tau * tau) / (h * h) * (solution[i-1][j+1] - 2 * solution[i-1][j] + solution[i-1][j-1]) + 2 * solution[i-1][j] - solution[i-2][j];
		}
	}
	return solution;
}
#endif //WAVE_EQ_H
