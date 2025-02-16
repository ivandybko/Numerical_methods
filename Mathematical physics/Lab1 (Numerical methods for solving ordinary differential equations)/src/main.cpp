#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
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
		std::cerr << "Error opening file!" << std::endl;
		return;
	}

	size_t rows = data.size();
	size_t cols = data[0].size();

	for (size_t j = 0; j < cols; ++j) {
		for (size_t i = 0; i < rows; ++i) {
			outFile << data[i][j] << " ";
		}
		outFile << "\n";
	}

	outFile.close();
}

template <typename T>
void AitkensProcess(T q, T tests, std::function<std::vector<std::vector<T>>(const std::vector<std::function<T(T, const std::vector<T>&)>>&, T, const std::vector<T>&, T, T)> solver, const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, T t_end, T tau) {
	size_t n = f.size();
	std::vector<std::vector<std::vector<T>>> data;
	// Для удобства сохраним используемые значения tau для каждой итерации
	std::vector<T> tau_vals;

	// Получаем решения для разных tau: tau_i = tau * pow(q, i)
	for (size_t i = 0; i < tests; ++i) {
		T current_tau = tau * std::pow(q, i);
		tau_vals.push_back(current_tau);
		std::vector<std::vector<T>> result = solver(f, t0, u0, t_end, current_tau);
		data.push_back(std::move(result));
	}
	// for (size_t i = 0; i < tests; ++i) {
	// 	T current_tau = tau * pow(q, i);
	// 	// std::cout << "Шаг tau_" << i << " = " << current_tau << '\n';
	// 	std::vector<std::vector<T>> result = solver(f, t0, u0, t_end, current_tau);
	// 	data.push_back(std::move(result));
	// }
	// std::vector<std::vector<T>> errors(tests);
	// for (size_t i = 0; i < tests; ++i) {
	// 	// Ошибку вычисляем для каждого компонента решения
	// 	errors[i].resize(n);
	// 	for (size_t j = 0; j < n; ++j) {
	// 		errors[i][j] = std::abs(data[i][j].back() - data.back()[j].back());
	// 	}
	// }
	T t_mid = (t0 + t_end) / 2;

	// Вычисляем ошибки для каждой переменной и для каждого набора tau.
	// Ошибка считается как разность между значением, полученным при данном tau,
	// и значением, полученным на самой точной (финальной) сетке, в точке t_mid.
	std::vector<std::vector<T>> errors(tests);
	for (size_t i = 0; i < tests; ++i) {
		// Для i-го набора tau вычисляем соответствующий индекс.
		// При равномерной сетке: t = t0 + index * (tau * pow(q, i))
		size_t idx_i = static_cast<size_t>(std::round((t_mid - t0) / tau_vals[i]));

		// Убедимся, что индекс не выходит за границы (если решатель корректно подбирает шаг, проблем не будет)
		errors[i].resize(n);
		for (size_t j = 0; j < n; ++j) {
			// Для решения с наименьшим tau (наиболее точного) вычисляем соответствующий индекс
			size_t idx_fine = static_cast<size_t>(std::round((t_mid - t0) / tau_vals.back()));
			errors[i][j] = std::abs(data[i][j][idx_i] - data.back()[j][idx_fine]);
		}
	}

	for (size_t i = 0; i < tests - 2; ++i) {
		std::cout << "Порядок сходимости p_" << i << ": ";
		for (size_t j = 0; j < errors[i].size(); ++j) {
			T p = log(errors[i + 1][j] / errors[i][j]) / log(q);
			std::cout << "[переменная " << j << ": " << p << "] ";
		}
		std::cout << '\n';
	}
}
template <typename T>

void AitkensProcess(T q, T tests, std::function<std::vector<std::vector<T>>(const std::vector<std::function<T(T, const std::vector<T>&)>>&, T, const std::vector<T>&, T, T, bool)> solver, const std::vector<std::function<T(T, const std::vector<T>&)>>& f,
	T t0, const std::vector<T>& u0, T t_end, T tau, bool predictor) {
	size_t n = f.size();
	std::vector<std::vector<std::vector<T>>> data;
	std::vector<T> tau_vals;
	for (size_t i = 0; i < tests; ++i) {
		T current_tau = tau * std::pow(q, i);
		tau_vals.push_back(current_tau);
		std::vector<std::vector<T>> result = solver(f, t0, u0, t_end, current_tau,predictor);
		data.push_back(std::move(result));
	}
	T t_mid = (t0 + t_end) / 2;

	std::vector<std::vector<T>> errors(tests);
	for (size_t i = 0; i < tests; ++i) {
		size_t idx_i = static_cast<size_t>(std::round((t_mid - t0) / tau_vals[i]));
		errors[i].resize(n);
		for (size_t j = 0; j < n; ++j) {
			size_t idx_fine = static_cast<size_t>(std::round((t_mid - t0) / tau_vals.back()));
			errors[i][j] = std::abs(data[i][j][idx_i] - data.back()[j][idx_fine]);
		}
	}
	for (size_t i = 0; i < tests - 2; ++i) {
		std::cout << "Порядок сходимости p_" << i << ": ";
		for (size_t j = 0; j < errors[i].size(); ++j) {
			T p = log(errors[i + 1][j] / errors[i][j]) / log(q);
			std::cout << "[переменная " << j << ": " << p << "] ";
		}
		std::cout << '\n';
	}
}

int main()
{
	// std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
	// 	[](double t, const std::vector<double>& u) { return 10*(u[1]-u[0]);},
	// 	[](double t, const std::vector<double>& u) { return u[0]*(28-u[2]) - u[1] ;},
	// 	[](double t, const std::vector<double>& u) { return  u[0]*u[1]-8.0/3.0*u[2]; }
	// };
	// std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
	// 	[](double t, const std::vector<double>& u) { return u[1];},
	// 	[](double t, const std::vector<double>& u) { return -20.0/0.3*u[0] ; }
	// };
	// std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
	// 	[](double t, const std::vector<double>& u) { return 2*u[0]+u[1]*u[1] -1;},
	// 	[](double t, const std::vector<double>& u) { return 6* u[0] -u[1]*u[1] +1 ; }
	// };
	std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
		[](double t, const std::vector<double>& u) { return 1-u[0]*u[0]-u[1]*u[1];},
		[](double t, const std::vector<double>& u) { return 2*u[0] ; }
	};

	std::vector<double> u0 = {0.0, 0.0}; double tend=1;
	std::cout << "Runge-Kutta 4th order method" << std::endl;
	std::vector<std::vector<double>> result = runge_kutta4(funcs, 0.0, u0, tend, 0.01);
	exportData(result, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/runge.txt");
	for (const auto & i : result)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Explicit Euler method" << std::endl;
	std::vector<std::vector<double>> result2 = explicit_euler(funcs, 0.0, u0, tend, 0.01);
	exportData(result2, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/explicit_euler.txt");
	for (const auto & i : result2)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Implicit Euler method" << std::endl;
	std::vector<std::vector<double>> result3 = implicit_euler(funcs, 0.0, u0, tend, 0.01);
	exportData(result3, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/implicit_euler.txt");
	for (const auto & i : result3)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Trapezoidal rule method" << std::endl;
	std::vector<std::vector<double>> result4 = trapezoidal_rule_method(funcs, 0.0, u0, tend, 0.01);
	exportData(result4, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/symmetric.txt");
	for (const auto & i : result4)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Adams-Bashforth 4th order method" << std::endl;
	std::vector<std::vector<double>> result5 = adams_bashforth(funcs, 0.0, u0, tend, 0.01);
	exportData(result5, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/adams_bashforth.txt");
	for (const auto & i : result5)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Adams-Bashforth 4th order method with predictor-corrector" << std::endl;
	std::vector<std::vector<double>> result6 = adams_bashforth(funcs, 0.0, u0, tend, 0.01, true);
	exportData(result6, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test2/adams_bashforth_with_predictor_corrector.txt");
	for (const auto & i : result6)
	{
		std::cout << i << std::endl;
	}
	std::cout << "Явный метод Эйлера" << std::endl;
	AitkensProcess<double>(0.5,9,explicit_euler<double>,funcs, 0.0, u0, tend, 0.05);
	std::cout << "Неявный метод Эйлера" << std::endl;
	AitkensProcess<double>(0.5,9,implicit_euler<double>,funcs, 0.0, u0, tend, 0.05);
	std::cout << "Симметричная схема" << std::endl;
	AitkensProcess<double>(0.5,9,trapezoidal_rule_method<double>,funcs, 0.0, u0, tend, 0.05);
	std::cout << "Метод Рунге-Кутты" << std::endl;
	AitkensProcess<double>(0.5,9,runge_kutta4<double>,funcs, 0.0, u0, tend, 0.05);
	std::cout << "Метод Адамса" << std::endl;
	AitkensProcess<double>(0.5,9,adams_bashforth<double>,funcs, 0.0, u0, tend, 0.05, false);
	std::cout << "Метод предиктор-корректор" << std::endl;
	AitkensProcess<double>(0.5,9,adams_bashforth<double>,funcs, 0.0, u0, tend, 0.05, true);
}