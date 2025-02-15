#ifndef EULER_H
#define EULER_H
#include <functional>
#include <vector>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/Gauss.h"
#include "../../../Linear algebra/Lab5 (Solving nonlinear equations)/src/FindRoots.h"
template <typename T>
std::vector<std::vector<T>> explicit_euler(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, double t_end, double tau){
	size_t n = f.size();
	size_t N =std::round((t_end-t0)/tau+1);
	std::vector<std::vector<T>> solution(n,std::vector<T>(N));
	T t;
	for (size_t j = 0; j < n; j++) {solution[j][0]=u0[j];}
	for(int i = 0; i < N-1; i++)
	{
		t=t0+i*tau;
		std::vector<T> u_n(n);
		for(int j = 0; j < n; j++)
		{
			u_n[j]=solution[j][i];
		}
		for(int j = 0; j < n; j++)
		{
			solution[j][i+1]=solution[j][i]+tau*f[j](t,u_n);
		}
	}
	return solution;
}

template <typename T>
std::vector<std::vector<T>> implicit_euler(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, double t_end, double tau){
	size_t n = f.size();
	size_t N =std::round((t_end-t0)/tau+1);
	std::vector<std::vector<T>> solution(n,std::vector<T>(1)); T t;
	for (size_t j = 0; j < n; j++) {solution[j][0]=u0[j];}
	std::vector<T> y = u0;
	auto F = [&](const std::vector<T>& y_next) -> std::vector<T> {
		std::vector<T> residual(y_next.size());
		for (size_t i = 0; i < y_next.size(); ++i) {
			T f_val = f[i](t + tau, y_next);
			residual[i] = y_next[i] - y[i] - tau * f_val;
		}
		return residual;
	};
	for(int i = 0; i < N-1; i++)
	{
		t=t0+i*tau;
		std::vector<size_t> subdivisions(n);
		std::vector<std::pair<T, T>> bounds(n);
		for(int j = 0; j < n; j++)
		{
			bounds[j]=std::pair<T,T>(solution[j][i]-tau*y[j],solution[j][i]+tau*y[j]);
			subdivisions[j]=10;
		}
		while (true){
			auto data=newtonMethod<double>(F,tau,1000, Gauss<double>, bounds, subdivisions);
			if (data.empty())
			{
				for(int j = 0; j < n; j++)
				{
					bounds[j]=std::pair<T,T>(bounds[j].first-2*tau*solution[j][i],bounds[j].second+2*tau*solution[j][i]);
					subdivisions[j]*=2;
				}
			}
			else
			{
				for(int j = 0; j < n; j++)
				{
					solution[j].emplace_back(data[0][j]);
				}
				break;
			}
		}
	}
	return solution;
}

template <typename T>
std::vector<std::vector<T>> trapezoidal_rule_method(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, double t_end, double tau){
	size_t n = f.size();
	size_t N =std::round((t_end-t0)/tau+1);
	std::vector<std::vector<T>> solution(n,std::vector<T>(1)); T t;
	for (size_t j = 0; j < n; j++) {solution[j][0]=u0[j];}
	std::vector<T> y = u0;
	auto F = [&](const std::vector<T>& y_next) -> std::vector<T> {
		std::vector<T> residual(y_next.size());
		for (size_t i = 0; i < y_next.size(); ++i) {
			T f_val = f[i](t + tau, y_next);
			T f_val_cur = f[i](t, y);
			residual[i] = y_next[i] - y[i] - tau/2 * f_val - tau/2 * f_val_cur;
		}
		return residual;
	};
	for(int i = 0; i < N-1; i++)
	{
		t=t0+i*tau;
		std::vector<size_t> subdivisions(n);
		std::vector<std::pair<T, T>> bounds(n);
		for(int j = 0; j < n; j++)
		{
			bounds[j]=std::pair<T,T>(solution[j][i]-tau*y[j],solution[j][i]+tau*y[j]);
			subdivisions[j]=10;
		}
		while (true){
			auto data=newtonMethod<double>(F,tau,1000,Gauss<double>, bounds, subdivisions);
			if (data.empty())
			{
				for(int j = 0; j < n; j++)
				{
					bounds[j]=std::pair<T,T>(bounds[j].first-2*tau*solution[j][i],bounds[j].second+2*tau*solution[j][i]);
					subdivisions[j]*=2;
				}
			}
			else
			{
				for(int j = 0; j < n; j++)
				{
					solution[j].emplace_back(data[0][j]);
				}
				break;
			}
		}
	}
	return solution;
}

#endif //EULER_H
