#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab3 (Solving interpolation problems)/src/grid_generators.h"
#include "runge_kutta_4th_order.h"
#include "euler.h"
#include "adams_bashforth_4th_order.h"

void exportData(const std::vector<std::vector<double>>& data, const std::string& filename) {
	if (data.empty()) return;

	std::ofstream outFile(filename);
	if (!outFile) {
		std::cerr << "Error opening file!" << '\n';
		return;
	}

	size_t const rows = data.size();
	size_t const cols = data[0].size();

	for (size_t j = 0; j < cols; ++j) {
		for (size_t i = 0; i < rows; ++i) {
			outFile << std::setprecision(16) << data[i][j] << " ";
		}
		outFile << "\n";
	}

	outFile.close();
}

template <typename T>
void AitkensProcess(T q, T tests, std::function<std::vector<std::vector<T>>(const std::vector<std::function<T(T, const std::vector<T>&)>>&, T, const std::vector<T>&, T, T)> solver, const std::vector<std::function<T(T, const std::vector<T>&)>>& func,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau) {
	size_t equations_number = func.size();
	std::vector<std::vector<std::vector<T>>> data;
	std::vector<T> tau_vals;
	for (size_t i = 0; i < tests; ++i) {
		T current_tau = tau * std::pow(q, i);
		tau_vals.push_back(current_tau);
		std::vector<std::vector<T>> result = solver(func, t_start, initial_conditions, t_end, current_tau);
		data.push_back(std::move(result));
	}
	T t_mid = (t_start + t_end) / 2;
	std::vector<std::vector<T>> errors(tests);
	for (size_t i = 0; i < tests; ++i) {
		auto idx_i = static_cast<size_t>(std::round((t_mid - t_start) / tau_vals[i]));
		errors[i].resize(equations_number);
		for (size_t j = 0; j < equations_number; ++j) {
			auto idx_fine = static_cast<size_t>(std::round((t_mid - t_start) / tau_vals.back()));
			errors[i][j] = std::abs(data[i][j][idx_i] - data.back()[j][idx_fine]);
		}
	}

	for (size_t i = 0; i < tests - 2; ++i) {
		std::cout << "Order of convergence p_" << i << ": ";
		std::cout << log(compute2Norm(errors[i + 1]) / compute2Norm(errors[i])) / log(q);
		std::cout << '\n';
	}
}
template <typename T>

void AitkensProcess(T q, T tests, std::function<std::vector<std::vector<T>>(const std::vector<std::function<T(T, const std::vector<T>&)>>&, T, const std::vector<T>&, T, T, bool)> solver, const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t_start, const std::vector<T>& initial_conditions, T t_end, T tau, bool predictor) {
	size_t equations_number = f.size();
	std::vector<std::vector<std::vector<T>>> data;
	std::vector<T> tau_vals;
	for (size_t i = 0; i < tests; ++i) {
		T current_tau = tau * std::pow(q, i);
		tau_vals.push_back(current_tau);
		std::vector<std::vector<T>> result = solver(f, t_start, initial_conditions, t_end, current_tau,predictor);
		data.push_back(std::move(result));
	}
	T t_mid = (t_start + t_end) / 2;

	std::vector<std::vector<T>> errors(tests);
	for (size_t i = 0; i < tests; ++i) {
		auto idx_i = static_cast<size_t>(std::round((t_mid - t_start) / tau_vals[i]));
		errors[i].resize(equations_number);
		for (size_t j = 0; j < equations_number; ++j) {
			auto idx_fine = static_cast<size_t>(std::round((t_mid - t_start) / tau_vals.back()));
			errors[i][j] = std::abs(data[i][j][idx_i] - data.back()[j][idx_fine]);
		}
	}
	for (size_t i = 0; i < tests - 2; ++i) {
		std::cout << "Order of convergence p_" << i << ": ";
		std::cout << log(compute2Norm(errors[i + 1]) / compute2Norm(errors[i])) / log(q);
		std::cout << '\n';
	}
}

int main()
{
	//Examples of systems
	std::vector<std::function<double(double, const std::vector<double>&)>> const pendulum = {
		[](double t, const std::vector<double>& u) { return u[1];},
		[](double t, const std::vector<double>& u) { return -20.0/0.3*u[0] ; }
	};
	std::vector<std::function<double(double, const std::vector<double>&)>> const test1 = {
		[](double t, const std::vector<double>& u) { return 2*u[0]+u[1]*u[1] -1;},
		[](double t, const std::vector<double>& u) { return 6* u[0] -u[1]*u[1] +1 ; }
	};
	std::vector<std::function<double(double, const std::vector<double>&)>> const test2 = {
		[](double t, const std::vector<double>& u) { return 1-u[0]*u[0]-u[1]*u[1];},
		[](double t, const std::vector<double>& u) { return 2*u[0] ; }
	};
	std::vector<std::function<double(double, const std::vector<double>&)>> const test3 = {
		[](double t, const std::vector<double>& u) { return 10*(u[1]-u[0]);},
		[](double t, const std::vector<double>& u) { return u[0]*(28-u[2]) - u[1] ;},
		[](double t, const std::vector<double>& u) { return  u[0]*u[1]-8.0/3.0*u[2]; }
	};
	std::vector<std::function<double(double, const std::vector<double>&)>> const bautin_model = {
		[](double t, const std::vector<double>& u) {
			double const a = 0.4 + 0.3 * std::sin(0.5 * t);
			double const E1 = 0.1;
			return a * u[1] + E1 * u[0] - std::pow(u[0], 3) - u[0] * std::pow(u[1], 2);
		},
		[](double t, const std::vector<double>& u) {
			double const a = 0.4 + 0.3 * std::sin(0.5 * t);
			double const E2 = 0.2;
			return -a * u[0] + E2 * u[1] - std::pow(u[1], 3) - u[1] * std::pow(u[0], 2);
		}
	};

	std::vector<std::function<double(double, const std::vector<double>&)>> const lotka_volterra_model = {
		[](double t, const std::vector<double>& u) {
			double const r1 = 0.4;
			double const b11 = 0.05;
			double const b12 = 0.1;
			return r1 * u[0] - b11 * u[0] * u[0] - b12 * u[0] * u[1];
		},
		[](double t, const std::vector<double>& u) {
			double const r2 = 0.1;
			double const b21 = 0.08;
			double const b22 = 0.003;
			return -r2 * u[1] - b22 * u[1] * u[1] + b21 * u[0] * u[1];
		}
	};
	//Setting initial conditions and end point
	std::vector<double> const u0 = {0.0,0.0}; double const tend=4;

	std::cout << "Explicit Euler method" << std::endl;
	std::vector<std::vector<double>> result2 = explicit_euler(test2, 0.0, u0, tend, 0.05);
	// exportData(result2, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/explicit_euler.txt");
	for (const auto & i : result2)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Implicit Euler method" << std::endl;
	std::vector<std::vector<double>> result3 = implicit_euler(test2, 0.0, u0, tend, 0.05);
	// exportData(result3, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/implicit_euler.txt");
	for (const auto & i : result3)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Trapezoidal rule method" << std::endl;
	std::vector<std::vector<double>> result4 = trapezoidal_rule_method(test2, 0.0, u0, tend, 0.05);
	// exportData(result4, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/symmetric.txt");
	for (const auto & i : result4)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Runge-Kutta 4th order method" << std::endl;
	std::vector<std::vector<double>> result = runge_kutta4(test2, 0.0, u0, tend, 0.1, false);
	// exportData(result, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/runge.txt");
	for (const auto & i : result)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Runge-Kutta 4th order method with step" << '\n';
	std::vector<std::vector<double>> const result_step = runge_kutta4(test2, 0.0, u0, tend, 0.1,true);
	// exportData(result_step, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/runge_step.txt");
	for (const auto & i : result_step)
	{
		std::cout << i << '\n';
	}
	std::cout << "Adams-Bashforth 4th order method" << std::endl;
	std::vector<std::vector<double>> result5 = adams_bashforth(test2, 0.0, u0, tend, 0.05);
	// exportData(result5, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/adams_bashforth.txt");
	for (const auto & i : result5)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Adams-Bashforth 4th order method with predictor-corrector" << std::endl;
	std::vector<std::vector<double>> result6 = adams_bashforth(test2, 0.0, u0, tend, 0.05, true);
	// exportData(result6, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/adams_bashforth_with_predictor_corrector.txt");
	for (const auto & i : result6)
	{
		std::cout << i << std::endl;
	}


	std::cout << "Explicit Euler method" << std::endl;
	AitkensProcess<double>(0.5,9,explicit_euler<double>, test2, 0.0, u0, tend, 0.05);
	std::cout << "Implicit Euler method" << std::endl;
	AitkensProcess<double>(0.5,9,implicit_euler<double>,test2, 0.0, u0, tend, 0.05);
	std::cout << "Trapezoidal rule method" << std::endl;
	AitkensProcess<double>(0.5,9,trapezoidal_rule_method<double>,test2, 0.0, u0, tend, 0.05);
	std::cout << "Runge-Kutta 4th order method" << std::endl;
	AitkensProcess<double>(0.5,9,runge_kutta4<double>,test2, 0.0, u0, tend, 0.05, false);
	std::cout << "Adams-Bashforth 4th order method" << std::endl;
	AitkensProcess<double>(0.5,9,adams_bashforth<double>,test2, 0.0, u0, tend, 0.05, false);
	std::cout << "Adams-Bashforth 4th order method with predictor-corrector" << std::endl;
	AitkensProcess<double>(0.5,9,adams_bashforth<double>,test2, 0.0, u0, tend, 0.05, true);
}