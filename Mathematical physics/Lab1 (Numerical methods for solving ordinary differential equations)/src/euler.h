#ifndef EULER_H
#define EULER_H
#include <functional>
#include <vector>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/Gauss.h"
#include "../../../Linear algebra/Lab5 (Solving nonlinear equations)/src/FindRoots.h"
template <typename T>
std::vector<std::vector<T>> explicit_euler(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& func,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau){
	size_t equations_number = func.size();
	size_t points_number =std::round((t_end-t_start)/tau+1);
	std::vector<std::vector<T>> solution(equations_number,std::vector<T>(points_number));
	T t;
	for (size_t j = 0; j < equations_number; j++) {solution[j][0]=initial_conditions[j];}
	for(int i = 0; i < points_number-1; i++)
	{
		t=t_start+i*tau;
		std::vector<T> u_n(equations_number);
		for(int j = 0; j < equations_number; j++)
		{
			u_n[j]=solution[j][i];
		}
		for(int j = 0; j < equations_number; j++)
		{
			solution[j][i+1]=solution[j][i]+tau*func[j](t,u_n);
		}
	}
	return solution;
}

template <typename T>
std::vector<std::vector<T>> implicit_euler(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& func,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau){
	size_t equations_number = func.size();
	size_t points_number =std::round((t_end-t_start)/tau+1);
	const int max_iterations = 1000;
	const int boundary_expansion = 10;
	std::vector<std::vector<T>> solution(equations_number,std::vector<T>(1)); T t;
	for (size_t j = 0; j < equations_number; j++) {solution[j][0]=initial_conditions[j];}
	std::vector<T> y = initial_conditions;
	auto residual_function = [&](const std::vector<T>& y_next) -> std::vector<T> {
		std::vector<T> residual(y_next.size());
		for (size_t i = 0; i < y_next.size(); ++i) {
			T f_val = func[i](t + tau, y_next);
			residual[i] = y_next[i] - y[i] - tau * f_val;
		}
		return residual;
	};
	for(int i = 0; i < points_number-1; i++)
	{
		t=t_start+i*tau;
		std::vector<size_t> subdivisions(equations_number);
		std::vector<std::pair<T, T>> bounds(equations_number);
		for(int j = 0; j < equations_number; j++)
		{
			bounds[j]=std::pair<T,T>(solution[j][i]-boundary_expansion*tau*(abs(y[j]+1)),solution[j][i]+boundary_expansion*tau*(abs(y[j]+1)));
			subdivisions[j]=boundary_expansion;
		}
		while (true){
			auto data=newtonMethod<T>(residual_function,tau*tau*tau, max_iterations, Gauss<double>, bounds, subdivisions);
			if (data.empty())
			{
				for(int j = 0; j < equations_number; j++)
				{
					bounds[j]=std::pair<T,T>(bounds[j].first-boundary_expansion*tau*(abs(solution[j][i])+1),bounds[j].second+boundary_expansion*tau*(abs(solution[j][i])+1));
					subdivisions[j]*=2;
				}
			}
			else
			{
				for(int j = 0; j < equations_number; j++)
				{
					solution[j].emplace_back(data[0][j]);
					y[j]=solution[j][i+1];
				}
				break;
			}
		}
	}
	return solution;
}

template <typename T>
std::vector<std::vector<T>> trapezoidal_rule_method(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& func,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau){
	const int max_iterations = 1000;
	const int boundary_expansion = 10;
	int const grid_refinement = 10;
	size_t equations_number = func.size();
	size_t points_number =std::round((t_end-t_start)/tau+1);
	std::vector<std::vector<T>> solution(equations_number,std::vector<T>(1)); T t;
	for (size_t j = 0; j < equations_number; j++) {solution[j][0]=initial_conditions[j];}
	std::vector<T> y = initial_conditions;
	auto residual_function = [&](const std::vector<T>& y_next) -> std::vector<T> {
		std::vector<T> residual(y_next.size());
		for (size_t i = 0; i < y_next.size(); ++i) {
			T f_val = func[i](t + tau, y_next);
			T f_val_cur = func[i](t, y);
			residual[i] = y_next[i] - y[i] - tau/2 * f_val - tau/2 * f_val_cur;
		}
		return residual;
	};
	for(int i = 0; i < points_number-1; i++)
	{
		t=t_start+i*tau;
		std::vector<size_t> subdivisions(equations_number);
		std::vector<std::pair<T, T>> bounds(equations_number);
		for(int j = 0; j < equations_number; j++)
		{
			bounds[j]=std::pair<T,T>(solution[j][i]-boundary_expansion*tau*(abs(y[j]+1)),solution[j][i]+boundary_expansion*tau*(abs(y[j]+1)));
			subdivisions[j]=grid_refinement;
		}
		while (true){
			auto data=newtonMethod<T>(residual_function,tau*tau*tau,max_iterations,Gauss<double>, bounds, subdivisions);
			if (data.empty())
			{
				for(int j = 0; j < equations_number; j++)
				{
					bounds[j]=std::pair<T,T>(bounds[j].first-boundary_expansion*tau*(abs(solution[j][i])+1),bounds[j].second+boundary_expansion*tau*(abs(solution[j][i])+1));
					subdivisions[j]*=2;
				}
			}
			else
			{
				for(int j = 0; j < equations_number; j++)
				{
					solution[j].emplace_back(data[0][j]);
					y[j]=solution[j][i+1];
				}
				break;
			}
		}
	}
	return solution;
}

#endif //EULER_H
