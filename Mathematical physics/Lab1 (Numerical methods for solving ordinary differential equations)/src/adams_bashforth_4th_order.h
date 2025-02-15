#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H
#include <functional>
#include <vector>
#include "runge_kutta_4th_order.h"
template<typename T>
std::vector<std::vector<T>> adams_bashforth(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, T t_end, T tau, bool predictor_corrector = false){
	size_t n = f.size(); T t;
	size_t N =(t_end-t0)/tau+1;
	auto solution = runge_kutta4(f, t0, u0, t0+3*tau, tau);
	for (auto &sol : solution){sol.resize(N);}
	for(int i = 3; i < N-1; i++)
	{
		t=t0+i*tau;
		std::vector<T> u_n(n), u_n1(n), u_n2(n), u_n3(n);
		for(int j = 0; j < n; j++)
		{
			u_n[j]=solution[j][i];
			u_n1[j]=solution[j][i-1];
			u_n2[j]=solution[j][i-2];
			u_n3[j]=solution[j][i-3];

		}
		for(int j = 0; j < n; j++)
		{
			solution[j][i + 1] = solution[j][i] + tau/24.0*(55.0* f[j](t, u_n)-59.0*f[j](t - tau, u_n1)
			+37.0*f[j](t - 2*tau, u_n2)-9.0*f[j](t - 3* tau, u_n3));
		}
		if ( predictor_corrector )
		{
			std::vector<T> u_cor(n);
			for(int j = 0; j < n; j++)
			{
				u_cor[j]=solution[j][i+1];
			}
			for(int j = 0; j < n; j++)
			{
				solution[j][i + 1] = solution[j][i] + tau/24.0*(9.0* f[j](t+tau, u_cor)+19.0*f[j](t, u_n)
				-5.0*f[j](t - tau, u_n1)+f[j](t - 2* tau, u_n2));
			}
		}

	}
	return solution;
}
#endif //ADAMS_BASHFORTH_H
