#include <iostream>
#include <fstream>
#include "heat_1d_eq_spatial_k.h"
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
			outFile << std::setprecision(20) << grid[k].first << " " << grid[k].second << " " << data[j][i] << '\n';
			k++;
		}
		k=k+pass*cols;
		j+=pass;
	}

	outFile.close();
}

double u(double x, double t, double L) {
	const double pi_ln2 = M_PI / std::log(2);
	double factor = 1.0 + x / L;
	return std::pow(factor, -0.5) * std::sin(pi_ln2 * std::log(factor)) *
		   std::exp(-t / (L * L) * (pi_ln2 * pi_ln2 + 0.25));
}

std::vector<double> error(double h, double tau, double sigma, double L, double t0, std::function<double(std::vector<double>)>  u0, std::function<double(std::vector<double>)> K, std::function<double(std::vector<double>)> initial_cond)
{
	auto data = heat_eq_spatial<double>(1, 1, 1.0, 0.1, 0.5, 2.0/3.0, 10.0, t0, L, 0.0, K,initial_cond, u0,false,  u0, false, tau, h, 0.5);
 	size_t const rows = data.size();
	size_t const cols = data[0].size();
	int k = 0;
	std::vector<double> max{0,0,0,0,0};
	for (size_t j = 0; j < rows; ++j) {
		for (size_t i = 0; i < cols; ++i) {
			double diff = std::abs(data[j][i]-u(i*h, j*tau, L));
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
std::function<double(std::vector<double>)> u0 = [](std::vector<double> data){return 0.2;};
std::function<double(std::vector<double>)> u2 = [](std::vector<double> data){return 10.0;};
std::function<double(std::vector<double>)> const P = [](std::vector<double> data){return 0.0;};
std::function<double(std::vector<double>)> const P1 = [](std::vector<double> data) {
	double t = data[0];
	double Q = data[1];
	double t0 = data[2];
	if (0 <= t && t <= t0) {
		return Q;
	}
	return 0.0;
};
std::function<double(std::vector<double>)> const P4 = [](std::vector<double> data) {
	double t = data[0];
	double Q = data[1];
	double t0 = data[2];
	if (0 <= t && t <= 0.5 * t0) {
		return 2 * Q * t;
	}
	if (0.5 * t0 < t && t < t0) {
			return 2 * Q * (t0 - t);
	}
	return 0.0;
};
// std::function<double(std::vector<double>)> const P1 = [](std::vector<double> data) {
// 	return 1.0;
// };
std::function<double(std::vector<double>)> const initial_cond = [](std::vector<double> data) {
	double x = data[0];
	double l = data[1];
	double u0 = data[2];
	return u0 + x*(l-x);
};
std::function<double(std::vector<double>)> const initial_cond0 = [](std::vector<double> data) {
	// return data[0];
	return 0.2;
};
std::function<double(std::vector<double>)> const initial_cond_L = [](std::vector<double> data) {
	double x = data[0];
	double L = data[1];
	return std::pow(1.0 + x / L, -0.5) *
	   std::sin((M_PI / std::log(2.0)) * std::log(1.0 + x / L));
};
std::function<double(std::vector<double>)> const K_L = [](std::vector<double> data) {
	double x = data[0];
	double L = data[1];
	return std::pow(1.0 + x / L, 2);
};
std::function<double(std::vector<double>)> const K = [](std::vector<double> data) {
	double x = data[0];
	double x1 = data[1];
	double x2 = data[2];
	double k1 = data[3];
	double k2 = data[4];
	double L = data[5];
	if (0 <= x && x <= x1) {
		return k1;
	}
	if (x1 < x && x < x2) {
		return k1 * (x - x2) / (x1 - x2) + k2 * (x - x1) / (x2 - x1);
	}
	if (x2 <= x && x <= L) {
		return k2;
	}
	std::cerr << "Function K(x) out of range!" << std::endl;
	return 0.0;
};
std::function<double(std::vector<double>)> const K_nonlinear = [](std::vector<double> data) {
	double u = data[0];
	double alpha = data[1];
	double beta = data[2];
	double gamma = data[3];
	return alpha + beta*std::pow(u, gamma);
};

std::function<double(std::vector<double>)> const left_nonlinear = [](std::vector<double> data)
{
	double t = data[0];
	double beta = data[1];
	double gamma = data[2];
	double c = data[3];
	return pow(((gamma * c * c) / beta), 1.0 / gamma) * pow(t, 1.0 / gamma);
};
std::function<double(std::vector<double>)> const right_nonlinear = [](std::vector<double> data)
{
	return 0;
};
std::function<double(std::vector<double>)> const K1 = [](std::vector<double> data) {
	return 100.0;
};

int main()
{


	double tau=0.01;
	double h=0.01;
	double t0{0.5}, L{1.0};

	auto sol0 = heat_eq_spatial<double>(0.25, 2, 2, 0.5, 0.5, 2.0/3.0, 10.0, t0, L, 0.2, K,initial_cond0, initial_cond0,false,  P4, true, tau, h, 0.5);
	auto grid = generate_uniform_grid<double>(0, t0, t0/tau+1, 0, L, L/h+1);
	exportData(grid, sol0, 24, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab2 (Numerical solution of boundary value problems for the one-dimensional heat equation)/data/test21.txt");
	// auto sol5 = heat_eq_spatial<double>(1, 1, 1.0, 0.1, 0.5, 2.0/3.0, 10.0, t0, L, 0.0, K_L,initial_cond_L, u0,false,  u0, false, tau, h, 0.5);
	// exportData(grid, sol5, 10, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab2 (Numerical solution of boundary value problems for the one-dimensional heat equation)/data/test5.txt");
	// auto sol1 = heat_eq_spatial<double>(1, 1, 1.0, 0.1, 0.5, 2.0/3.0, 10.0, t0, L, 0.0, K_L,initial_cond_L, u0,false,  u0, false, tau, h, 1.0);
	// exportData(grid, sol1, 10, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab2 (Numerical solution of boundary value problems for the one-dimensional heat equation)/data/test1.txt");

	// auto sol8 = heat_eq_spatial<double>(1.25, 4, 1.5, 0.1, 0.4, 0.5, 10.0, t0, L, 0.1, K,u0, u0,false,  P4, true, tau, h, 1.0);
	// auto grid = generate_uniform_grid<double>(0, t0, t0/tau+1, 0, L, L/h+1);
	// exportData(grid, sol8, 10, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab2 (Numerical solution of boundary value problems for the one-dimensional heat equation)/data/test21.txt");
	for (auto &x : sol0)
	{
		std::cout << x << std::endl;
	}
	std::ofstream outFile("/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab2 (Numerical solution of boundary value problems for the one-dimensional heat equation)/data/err.txt");
	for (int i=2; i<4 ; i++)
	{
		auto err = error(h, tau, 0.5, L, t0, u0, K_L, initial_cond_L);
		for (auto elem : err)
		{
			outFile << elem << ' ';
		}
		outFile << '\n';
		std::cout << h << '\n';
		// h = 1.0/(i*2.0);
		h/=2;
		tau/=2;
	}
	return 0;
}
