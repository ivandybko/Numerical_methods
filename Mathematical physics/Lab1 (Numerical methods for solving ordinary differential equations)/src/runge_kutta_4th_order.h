
#ifndef RUNGE_KUTTA_4TH_ORDER_CPP_H
#define RUNGE_KUTTA_4TH_ORDER_CPP_H
#include <functional>
#include <vector>
#include <chrono>


template<typename T>
std::vector<std::vector<T>> runge_kutta4(
	const std::vector<std::function<T(T, const std::vector<T>&)>>& func,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau, bool RungeStep = true){
	size_t equations_number = func.size();
	T begin_tau = tau;
	std::ofstream outFile("/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/lotka_volterra/step.txt");
	// size_t N =(t_end-t0)/tau+1;
	T tolerance = 1e-6;
	std::vector<T> k1(equations_number);
	std::vector<T> k2(equations_number);
	std::vector<T> k3(equations_number);
	std::vector<T> k4(equations_number);
	T t=t_start;
	std::vector<std::vector<T>> solution(equations_number+1,std::vector<T>(1));
	solution[0][0] = t_start;
	for (size_t j = 0; j < equations_number; j++) {solution[j+1][0]=initial_conditions[j];}
	int i=0;
	auto start_time = std::chrono::steady_clock::now();
	auto last_print_time = start_time;
	while (t < t_end){
			t+=tau;
			solution[0].push_back(t);
			std::vector<T> u_n(equations_number);
			for(int j = 0; j < equations_number; j++)
			{
				u_n[j]=solution[j+1][i];
			}
			for(int j = 0; j < equations_number; j++)
			{
				k1[j]=func[j](t,u_n);
			}
			for(int j = 0; j < equations_number; j++)
			{
				k2[j]=func[j](t+0.5*tau,u_n+0.5*tau*k1);
			}
			for(int j = 0; j < equations_number; j++)
			{
				k3[j]=func[j](t+0.5*tau,u_n+0.5*tau*k2);
			}
			for(int j = 0; j < equations_number; j++)
			{
				k4[j]=func[j](t+tau,u_n+tau*k3);
			}
			for(int j = 0; j < equations_number; j++)
			{
				u_n[j] = solution[j+1][i]+tau/6.0*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
			}
			T diff_max{0};
			if(RungeStep)
			{
				std::vector<T> u_n_half(equations_number);
				T t_half=t-tau/2.0;
				T tau_half=tau/2.0;
				for(int j = 0; j < equations_number; j++)
				{
					u_n_half[j]=solution[j+1][i];
				}
				for(int j = 0; j < equations_number; j++)
				{
					k1[j]=func[j](t_half,u_n);
				}
				for(int j = 0; j < equations_number; j++)
				{
					k2[j]=func[j](t_half+0.5*tau_half,u_n+0.5*tau_half*k1);
				}
				for(int j = 0; j < equations_number; j++)
				{
					k3[j]=func[j](t_half+0.5*tau_half,u_n+0.5*tau_half*k2);
				}
				for(int j = 0; j < equations_number; j++)
				{
					k4[j]=func[j](t_half+tau_half,u_n+tau_half*k3);
				}
				for(int j = 0; j < equations_number; j++)
				{
					u_n_half[j] = solution[j+1][i]+tau/6.0*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
				}
				t_half+=tau;
				for(int j = 0; j < equations_number; j++)
				{
					k1[j]=func[j](t_half,u_n);
				}
				for(int j = 0; j < equations_number; j++)
				{
					k2[j]=func[j](t_half+0.5*tau_half,u_n+0.5*tau_half*k1);
				}
				for(int j = 0; j < equations_number; j++)
				{
					k3[j]=func[j](t_half+0.5*tau_half,u_n+0.5*tau_half*k2);
				}
				for(int j = 0; j < equations_number; j++)
				{
					k4[j]=func[j](t_half+tau_half,u_n+tau_half*k3);
				}
				for(int j = 0; j < equations_number; j++)
				{
					u_n_half[j] = solution[j+1][i]+tau/6.0*(k1[j]+2*k2[j]+2*k3[j]+k4[j]);
				}
				diff_max=compute2Norm(u_n-u_n_half)/15;
			}
			if (RungeStep and diff_max > 1.05*tolerance ){
				T tau_old = tau;
				t-=tau;
				tau=tau*std::pow(tolerance/diff_max,0.25);
				i--;
				solution[0].pop_back();
				u_n.resize(0);
			}
			if (RungeStep and compute2Norm(u_n) > 10000){std::cout << "Automatic step turned off" << '\n'; RungeStep = false; tau=begin_tau; }
			if (RungeStep and diff_max < tolerance/2)
			{
				T tau_old = tau;
				tau=tau*std::pow(tolerance/diff_max,0.25);
			}
			if (u_n.size()==equations_number)
			{
				for(int j = 0; j < equations_number; j++)
				{
					solution[j+1].push_back(u_n[j]);
				}
			}
			auto current_time = std::chrono::steady_clock::now();
			if (current_time - last_print_time >= std::chrono::minutes(1)) {
				std::cout << "Calculated " << (t / t_end * 100) << "%" << std::endl;
				last_print_time = current_time;
			}
			outFile << t << " " << tau << " " << diff_max << std::endl;
			i++;
	}
	outFile.close();
	return solution;
}

#endif //RUNGE_KUTTA_4TH_ORDER_CPP_H
