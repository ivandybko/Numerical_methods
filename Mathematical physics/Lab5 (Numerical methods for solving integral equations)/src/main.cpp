#include <iostream>
#include <fstream>
#include "quadrature_fredholm_solver.h"
#include "simple_iteration_fredholm_solver.h"
#include "degenerate_kernel_fredholm_solver.h"
#include "singular_fredholm_solver.h"
template <typename T>
void data_export(std::vector<T> solution, T a, T b, int grid_size, const std::string& filename)
{
	if (solution.empty()) return;
	std::ofstream outFile(filename);
	if (!outFile) {
		std::cerr << "Error opening file!" << '\n';
		return;
	}
	T h = (b - a) / (grid_size - 1);
	for (int i = 0; i < solution.size(); i++)
	{
		outFile << a+i*h << " " << solution[i] << '\n';
	}

	outFile.close();
}

template <typename T>
T error(std::function<T(T)> exact, std::vector<T> solution, T a, T b, int grid_size)
{
	T h = (b - a) / (grid_size - 1);
	T max = 0.0;
	for (int i = 0; i < solution.size(); i++)
	{
		if (std::abs(solution[i] - exact(a+i*h)) > max)
		{
			max = std::abs(solution[i] - exact(a+i*h));
		}
	}
	return max;
}

// Первое уравнение
std::function<double(double, double)> K_1 = [](double x, double s) {
	return 0.5 * (1 - x * std::cos(x * s));
};

std::function<double(double)> f_1 = [](double x) {
	return 0.5 * (1 + std::sin(x));
};

// Второе уравнение
std::function<double(double, double)> K_2 = [](double x, double s) {
	return 0.5 * (1 - x * std::cos(x * s));
};

std::function<double(double)> f_2 = [](double x) {
	return x * x + std::sqrt(x);
};

std::function<double(double, double)> K_3 = [](double x, double s) {
	return std::sin(x) * std::cos(s);
};

std::function<double(double)> f_3 = [](double x) {
	return x * x;
};

double kernel8(double x, double s) {
	return 1 - x * std::cos(x * s);
}

double f1_8(double x) {
	return std::sin(x);
}

double f2_8(double x) {
	return std::sin(x) - std::pow(x, 3);
}
// Первое уравнение
std::function<double(double, double)> K_ord = [](double x, double s) {
	return x*s;
};

std::function<double(double)> f_ord = [](double x) {
	return 2.0/3.0 * x;
};

std::function<double(double, double)> K_v = [](double x, double s) {
	return std::exp(x*s);
};

std::function<double(double)> f_v = [](double x) {
	return x;
};

std::function<double(double)> sol_v = [](double x) {
	return x - 2.0 / (std::exp(2) + 1) * std::exp(x);
};


std::function<double(double)> sol_ord = [](double x) {
	return x;
};

// phi_k(x) = x^k / k!
std::vector<std::function<double(double)>> phi = {
	[](double x) { return 1.0; },                                           // k = 0
	[](double x) { return x; },                                             // k = 1
	[](double x) { return -x*x / 2.0; },                                    // k = 2
	[](double x) { return -x*x*x / 6.0; },                                  // k = 3
	[](double x) { return x*x*x*x / 24.0; },                                // k = 4
	[](double x) { return x*x*x*x*x / 120.0; },                             // k = 5
	[](double x) { return -x*x*x*x*x*x / 720.0; },                          // k = 6
	[](double x) { return -x*x*x*x*x*x*x / 5040.0; },                       // k = 7
	[](double x) { return x*x*x*x*x*x*x*x / 40320.0; },                     // k = 8
	[](double x) { return x*x*x*x*x*x*x*x*x / 362880.0; },                 // k = 9
	[](double x) { return -x*x*x*x*x*x*x*x*x*x / 3628800.0; },             // k = 10
	[](double x) { return -x*x*x*x*x*x*x*x*x*x*x / 39916800.0; },          // k = 11
	[](double x) { return x*x*x*x*x*x*x*x*x*x*x*x / 479001600.0; }         // k = 12
};

// psi_k(s) = s^k
std::vector<std::function<double(double)>> psi = {
	[](double s) { return 1.0; },                                           // s^0
	[](double s) { return s; },
	[](double s) { return std::pow(s, 2); },
	[](double s) { return std::pow(s, 3); },
	[](double s) { return std::pow(s, 4); },
	[](double s) { return std::pow(s, 5); },
	[](double s) { return std::pow(s, 6); },
	[](double s) { return std::pow(s, 7); },
	[](double s) { return std::pow(s, 8); },
	[](double s) { return std::pow(s, 9); },
	[](double s) { return std::pow(s, 10); },
	[](double s) { return std::pow(s, 11); },
	[](double s) { return std::pow(s, 12); }
};

std::vector<std::function<double(double)>> phi_2 = {
	[](double x) { return 1.0; },                                           // k = 0

};

std::vector<std::function<double(double)>> psi_2 = {
	[](double s) { return 1.0; },                                           // s^0
};

std::function<double(double, double)> K_div = [](double x, double t) {
	return  std::cos(x * t) + std::sin(x*t);
};

std::function<double(double)> f_div = [](double x) -> double {
	// Проверка на недопустимые значения: x = -1 или x = 1 (деление на ноль)
	if (std::abs(x + 1.0) < 1e-12 || std::abs(x - 1.0) < 1e-12) {
		throw std::domain_error("Division by zero in function f(x): x = -1 or x = 1");
	}

	double term1 = -1.0 / (1.0 + x);
	double term2 = std::cos(x);
	double term3 = std::cos((1.0 + x) / 2.0) / (1.0 + x);
	double term4 = std::sin((1.0 - x) / 2.0) / (-1.0 + x);
	double term5 = std::sin(x);

	return term1 + term2 + term3 + term4 + term5;
};

std::function<double(double)> f_div_2 = [](double x) -> double {
	return x - 0.5;
};

std::function<double(double)> sol_div_e = [](double x) {
	return  std::cos(x) + std::sin(x);
};
// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main()
{
	int grid_size = 50;
	// auto sol = quadrature_fredholm_solve<double>(1, 0, 1, K_ord, f_ord, grid_size);
	// std::string path_1 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/ord.txt";
	// data_export<double>(sol, 0, 1, grid_size, path_1);
	// std::cout << error<double>(sol_ord, sol, 0, 1, grid_size);
	// auto sol2 = quadrature_fredholm_solve<double>(1, 0, 1, K_2, f_2, grid_size);
	// std::string path_2 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/2.txt";
	// data_export<double>(sol2, 0.1, 1, grid_size, path_2);

	// auto sol2_si = simple_iteration_fredholm_solve<double>(1, 0.0, 1, K_ord, f_ord, grid_size, 5, 1e-6);
	// std::cout << error<double>(sol_ord, sol2_si, 0, 1, grid_size);
	// std::string path_2_si = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/8_si.txt";
	// data_export<double>(sol2_si, 0, 1, grid_size, path_2_si);

	// auto sol_v_si = simple_iteration_fredholm_solve<double>(-1, 0.0, 1, K_v, f_v, grid_size, 10, 1e-6);
	// std::cout << error<double>(sol_v, sol_v_si, 0, 1, grid_size);
	// std::string path_v_si = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/v_si.txt";
	// data_export<double>(sol2_si, 0, 1, grid_size, path_2_si);

	// auto sol_div = degenerate_kernel_fredholm_solver<double>(1, 0.0 , 0.5, f_div, phi, psi, grid_size);
	// std::string path_div = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/div.txt";
	// data_export<double>(sol_div, 0, 0.5, grid_size, path_div);
	// std::cout << error<double>(sol_div_e, sol_div, 0, 0.5, grid_size) << std::endl;
	//
	//
	// auto sol_div_q = quadrature_fredholm_solve<double>(1, 0, 0.5, K_div, f_div, grid_size);
	// std::string path_div_q = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/div_q.txt";
	// data_export<double>(sol_div_q, 0, 0.5, grid_size, path_div_q);
	// std::cout << error<double>(sol_div_e, sol_div_q, 0, 0.5, grid_size) << std::endl; ;
	//
	//
	//
	// auto sol_div_2 = degenerate_kernel_fredholm_solver<double>(1, 0.0 , 1.0, f_div_2, phi_2, psi_2, grid_size);
	// std::string path_div_2 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/div_2.txt";
	// data_export<double>(sol_div_2, 0, 1.0, grid_size, path_div_2);
	// std::cout << error<double>(sol_ord, sol_div_2, 0, 1.0, grid_size) << std::endl;

	int N = 1000;
	std::function<double(double)> f_singular = [](double angle) {
		return std::cos(4*angle);
	};
	std::string filename = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/R.txt";
	std::ofstream outFile(filename);
	auto sol_singular = singular_fredholm_solver<double>(f_singular, N);
	for (int N = 1; N < 50; N++)
	{
		auto sol_singular = singular_fredholm_solver<double>(f_singular, N);
		outFile << N << " " << sol_singular.back() << " " << '\n';
	}
	outFile.close();
	// std::string path_singular = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/singular.txt";
	// data_export<double>(sol_singular, 0, 1, grid_size, path_singular);
	std::cout << sol_singular << '\n';
	// std::string order = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/order.txt";
	// std::ofstream outFile(order);
	// int iter = 1;
	// for (int i = 0; i < 20; i++)
	// {
	// 	auto sol2_si = simple_iteration_fredholm_solve<double>(1, 0.0, 1, K_ord, f_ord, grid_size, iter, 1e-6);
	// 	outFile << iter << " "<< error<double>(sol_ord, sol2_si, 0, 1, grid_size) << '\n';
	// 	iter++;
	// }
	// for (int i = 1; i < 12; i++)
	// {
	// 	std::vector<std::function<double(double)>> x(i+1), y(i+1);
	// 	for (int j = 0; j <= i; j++)
	// 	{
	// 		x[j] = phi[j];
	// 		y[j] = psi[j];
	// 	}
	// 	auto sol_div = degenerate_kernel_fredholm_solver<double>(1, 0.0 , 0.5, f_div, x, y, grid_size);
	// 	// std::string path_div = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/div.txt";
	// 	// data_export<double>(sol_div, 0, 0.5, grid_size, path_div);
	// 	std::cout << error<double>(sol_div_e, sol_div, 0, 0.5, grid_size) << std::endl;
	// }
	return 0;
}

// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.