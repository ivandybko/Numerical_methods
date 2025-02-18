
#ifndef RUNGE_KUTTA_4TH_ORDER_CPP_H
#define RUNGE_KUTTA_4TH_ORDER_CPP_H
#include <functional>
#include <vector>
#include <chrono>


template<typename T>
std::vector<std::vector<T>> runge_kutta4(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, T t_end, T tau, bool RungeStep = true){
	size_t n = f.size();
	T begin_tau = tau;
	// size_t N =(t_end-t0)/tau+1;
	T tolerance = 0.001;
	std::vector<T> k1(n),k2(n),k3(n),k4(n);
	T t=t0;
	std::vector<std::vector<T>> solution(n+1,std::vector<T>(1));
	solution[0][0] = t0;
	for (size_t j = 0; j < n; j++) {solution[j+1][0]=u0[j];}
	int i=0;
	auto start_time = std::chrono::steady_clock::now();
	auto last_print_time = start_time;
	while (t < t_end){
			t+=tau;
			solution[0].push_back(t);
			std::vector<T> u_n(n);
			for(int j = 0; j < n; j++)
			{
				u_n[j]=solution[j+1][i];
			}
			for(int j = 0; j < n; j++)
			{
				k1[j]=f[j](t,u_n);
			}
			for(int j = 0; j < n; j++)
			{
				k2[j]=f[j](t+0.5*tau,u_n+0.5*tau*k1);
			}
			for(int j = 0; j < n; j++)
			{
				k3[j]=f[j](t+0.5*tau,u_n+0.5*tau*k2);
			}
			for(int j = 0; j < n; j++)
			{
				k4[j]=f[j](t+tau,u_n+tau*k3);
			}
			std::vector<T> u, diff(n);
			for(int j = 0; j < n; j++)
			{
				auto test = solution[j+1][i]+tau/6.0*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
				diff[j] = std::abs((test - solution[j+1][i])/15);
				u.push_back(test);
			}
			T diff_max = *std::max_element(diff.begin(),diff.end());
			if (RungeStep and diff_max > 1.2*tolerance ){
				T tau_old = tau;
				t-=tau;
				tau=tau*std::pow(tolerance/diff_max,0.25);
				i--;
				solution[0].pop_back();
				u.resize(0);
			}
			if (RungeStep and compute2Norm(u) > 10000){std::cout << "Automatic step turned off" << '\n'; RungeStep = false; tau=begin_tau; }
			if (RungeStep and diff_max < tolerance/10)
			{
				T tau_old = tau;
			}
			if (u.size()==n)
			{
				for(int j = 0; j < n; j++)
				{
					solution[j+1].push_back(u[j]);
				}
			}
			auto current_time = std::chrono::steady_clock::now();
			if (current_time - last_print_time >= std::chrono::minutes(1)) {
				std::cout << "Calculated " << (t / t_end * 100) << "%" << std::endl;
				last_print_time = current_time;
			}
			i++;
	}
	return solution;
}

#endif //RUNGE_KUTTA_4TH_ORDER_CPP_H
