
#ifndef RUNGE_KUTTA_4TH_ORDER_CPP_H
#define RUNGE_KUTTA_4TH_ORDER_CPP_H
#include <functional>
#include <vector>
#include "../../../Linear algebra/Lab3 (Solving interpolation problems)/src/grid_generators.h"
template<typename T>
std::vector<std::vector<T>> RungeKutta4System(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, double t_end, double tau){
	size_t n = f.size();
	size_t N =(t_end-t0)/tau+1;
	std::vector<T> k1(n),k2(n),k3(n),k4(n),t(N);
	for (size_t i = 0; i < N; i++) {
		t[i]=t0+i*tau;
	}
	std::vector<std::vector<T>> solution(n,std::vector<T>(N));
	for (size_t j = 0; j < n; j++) {solution[j][0]=u0[j];}
	for(int i = 0; i < N-1; i++)
	{
		std::vector<T> u_n(n);
		for(int j = 0; j < n; j++)
		{
			u_n[j]=solution[j][i];
		}
		for(int j = 0; j < n; j++)
		{
			k1[j]=f[j](t[i],u_n);
		}
		for(int j = 0; j < n; j++)
		{
			k2[j]=f[j](t[i]+0.5*tau,u_n+0.5*tau*k1);
		}
		for(int j = 0; j < n; j++)
		{
			k3[j]=f[j](t[i]+0.5*tau,u_n+0.5*tau*k2);
		}
		for(int j = 0; j < n; j++)
		{
			k4[j]=f[j](t[i]+tau,u_n+tau*k3);
		}
		for(int j = 0; j < n; j++)
		{
			solution[j][i+1]=solution[j][i]+tau/6*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
		}
	}
	return solution;
}

#endif //RUNGE_KUTTA_4TH_ORDER_CPP_H
