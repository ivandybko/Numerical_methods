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
			outFile << std::setprecision(16) << grid[k].first << " " << grid[k].second << " " << data[j][i] << '\n';
			k++;
		}
		j+=pass;
	}

	outFile.close();
}

int main()
{
	std::function<double(double)> const u0 = [](double t){return 300;};
	std::function<double(double)> const u20 = [](double t){return 300;};
	std::function<double(double, double, double)> const P4 = [](double t, double Q, double t0) {
		if (0 <= t && t <= 0.5 * t0) {
			return 2 * Q * t;
		}
		if (0.5 * t0 < t && t < t0) {
				return 2 * Q * (t0 - t);
		}
		return 0.0;
	};
	std::function<double(double, double, double)> const initial_cond = [](double x, double u0, double l) {
		return u0 + x*(l-x);
	};
	std::function<double(double, double, double, double, double, double)> const K = [](double x, double x1, double x2, double k1, double k2, double L) {
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
	std::function<double(double, double, double, double, double, double)> const K1 = [](double x, double x1, double x2, double k1, double k2, double L) {
		return 5000;
	};
	double tau=10;
	double h=0.5;
	double t0{100000.0}, L{10.0};
	auto sol = heat_eq_spatial<double>(7800, 460, 1.0, 0.1, 0.5, 2.0/3.0, 10.0, t0, L, 300.0, K,initial_cond, u0, u0, tau, h, 1.0);
	auto grid = generate_uniform_grid<double>(0, t0, t0/tau+1, 0, L, L/h+1);
	exportData(grid, sol, 5, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab2 (Numerical solution of boundary value problems for the one-dimensional heat equation)/data/test.txt");
	for (auto &x : sol)
	{
		std::cout << x << std::endl;
	}
	return 0;
}
