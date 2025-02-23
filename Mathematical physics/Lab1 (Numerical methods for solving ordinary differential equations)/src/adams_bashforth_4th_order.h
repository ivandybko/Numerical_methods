#ifndef ADAMS_BASHFORTH_H
#define ADAMS_BASHFORTH_H
#include <functional>
#include <vector>
#include "runge_kutta_4th_order.h"
template<typename T>
std::vector<std::vector<T>> adams_bashforth(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& func,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau, bool predictor_corrector = false){
	size_t equations_number = func.size(); T t;
	size_t points_number =(t_end-t_start)/tau+1;
	auto solution = runge_kutta4(func, t_start, initial_conditions, t_start+3*tau, tau);
	for (auto &sol : solution){sol.resize(points_number);}
	for(int i = 3; i < points_number-1; i++)
	{
		t=t_start+i*tau;
		std::vector<T> u_n(equations_number);
		std::vector<T> u_n1(equations_number);
		std::vector<T> u_n2(equations_number);
		std::vector<T> u_n3(equations_number);
		for(int j = 0; j < equations_number; j++)
		{
			u_n[j]=solution[j][i];
			u_n1[j]=solution[j][i-1];
			u_n2[j]=solution[j][i-2];
			u_n3[j]=solution[j][i-3];

		}
		for(int j = 0; j < equations_number; j++)
		{
			solution[j][i + 1] = solution[j][i] + tau/24.0*(55.0* func[j](t, u_n)-59.0*func[j](t - tau, u_n1)
			+37.0*func[j](t - 2*tau, u_n2)-9.0*func[j](t - 3* tau, u_n3));
		}
		if ( predictor_corrector )
		{
			std::vector<T> u_cor(equations_number);
			for(int j = 0; j < equations_number; j++)
			{
				u_cor[j]=solution[j][i+1];
			}
			for(int j = 0; j < equations_number; j++)
			{
				solution[j][i + 1] = solution[j][i] + tau/24.0*(9.0* func[j](t+tau, u_cor)+19.0*func[j](t, u_n)
				-5.0*func[j](t - tau, u_n1)+func[j](t - 2* tau, u_n2));
			}
		}

	}
	return solution;
}
#endif //ADAMS_BASHFORTH_H
