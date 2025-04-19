#include <iostream>
#include <fstream>
#include "wave_eq.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab3 (Solving interpolation problems)/src/grid_generators.h"

void exportData(const std::vector<std::pair<double, double>>& grid, const std::vector<std::vector<double>>& data, int pass, const std::string& filename) {
	if (data.empty()) return;

	std::ofstream outFile(filename);
	if (!outFile) {
		std::cerr << "Error opening file!" << '\n';
		return;
	}

	size_t const rows = data.size();
	size_t const cols = data[0].size();
	int k = 0;
	for (size_t j = 0; j < rows; ++j) {
		for (size_t i = 0; i < cols; ++i) {
			outFile << std::setprecision(16) << grid[k].first << " " << grid[k].second << " " << data[j][i] << '\n';
			k++;
		}
		k=k+pass*cols;
		j+=pass;
	}

	outFile.close();
}

template <typename T>
std::vector<double> error(std::function<T(T,T)> sol, T tau, T h, T a, T t0, T left, T right, std::function<T(T)> initial_cond, std::function<T(T)> initial_speed, std::function<T(T)> boundary_left, std::function<T(T)> boundary_right, std::function<T(T)> initial_cond_d = nullptr)
{
	auto data = wave_eq(a, t0, left, right,initial_cond, initial_speed, boundary_left, boundary_right, tau, h, initial_cond_d);
	size_t const rows = data.size();
	size_t const cols = data[0].size();
	int k = 0;
	std::vector<double> max{0,0,0,0,0};
	for (size_t j = 0; j < rows; ++j) {
		for (size_t i = 0; i < cols; ++i) {
			double diff = std::abs(data[j][i]-sol(i*h, j*tau));
			if (diff > max[0])
			{
				max[0] = diff;
				max[1] = i*h;
				max[2] = j*tau;
				max[3] = h;
				max[4]= tau;
			}
			k++;
		}
	}
	return max;
}

std::function<double(double)> const initial_cond_ex1 = [](double x) {
	return std::sin(M_PI * x);
};

std::function<double(double)> const initial_cond_d_ex1 = [](double x) {
	return -M_PI*M_PI*std::sin(M_PI * x);
};

std::function<double(double)> const initial_cond_ex2 = [](double x) {
	return x*(1-x);
};

std::function<double(double)> const initial_speed_ex1 = [](double x) {
	return 0;
};
std::function<double(double)> const boundary_left_ex1 = [](double t) {
	return 0;
};
std::function<double(double)> const boundary_right_ex1 = [](double t) {
	return 0;
};

std::function<double(double)> const initial_cond_task_1 = [](double x) {
	if (x >= -1.0/3.0 && x <= 1.0/3.0)
	{
		return 1.0;
	}
	return 0.0;
};
std::function<double(double)> const initial_cond_task_2 = [](double x) {
	if (x >= -1.0/2.0 && x <= 1.0/2.0)
	{
		return 1 - 2* std::abs(x);
	}
	return 0.0;
};
double sol_ex1(double x, double t) {
	return std::sin(M_PI * x) * std::cos(M_PI * t);
}


int main()
{
	double tau=0.1;
	double h=0.1;
	double t0{10.0}, L{10.0};
	double left{0.0}, right{4.0 * M_PI};
	auto sol1 =  wave_eq<double>(1, t0, left, right, initial_speed_ex1, initial_speed_ex1, initial_cond_ex1, boundary_right_ex1, tau,h);
	auto grid = generate_uniform_grid<double>(0, t0, t0/tau+1, left, right, (right-left)/h+1);
	exportData(grid, sol1, 0, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab3 (Numerical solution of boundary value problems for the one-dimensional wave equation)/data/task3.txt");
	for (auto &x : sol1)
	{
		std::cout << x << std::endl;
	}
	// for (int i = 0; i < 5; ++i)
	// {
	// 	// auto sol1 =  wave_eq<double>(1, L, t0, initial_cond_ex1, initial_speed_ex1, boundary_left_ex1, boundary_right_ex1, tau,h);
	// 	// for (auto &x : sol1)
	// 	// {
	// 	// 	std::cout << x << std::endl;
	// 	// }
	// 	std::cout << "h=" <<  h << " tau=" << tau << '\n';
	// 	std::cout << error<double>(sol_ex1, tau, h, 1, t0, left, right, initial_cond_ex1, initial_speed_ex1, boundary_left_ex1, boundary_right_ex1, initial_cond_d_ex1) << '\n';
	// 	h/=2;
	// 	tau/=2;;
	// }
	return 0;
}
