#include <iostream>
#include <vector>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab3 (Solving interpolation problems)/src/grid_generators.h"
#include "runge_kutta_4th_order.h"
#include "euler.h"
#include <fstream>

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

int main()
{
	std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
		[](double t, const std::vector<double>& u) { return 2*u[0]+u[1]*u[1] -1;},
		[](double t, const std::vector<double>& u) { return 6* u[0] -u[1]*u[1] +1 ; }
	};
	// std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
	// 	[](double t, const std::vector<double>& u) { return -u[1];},
	// 	[](double t, const std::vector<double>& u) { return u[0] ; }
	// };

	std::vector<double> u0 = {1.0, 1.0}; double tend=0.5;
	std::vector<std::vector<double>> result = runge_kutta4(funcs, 0.0, u0, tend, 0.05);
	exportData(result, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test1/runge.txt");
	for (const auto & i : result)
	{
		std::cout << i << std::endl;
	}
	std::vector<std::vector<double>> result2 = explicit_euler(funcs, 0.0, u0, tend, 0.05);
	exportData(result, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test1/explicit_euler.txt");
	for (const auto & i : result2)
	{
		std::cout << i << std::endl;
	}
	std::vector<std::vector<double>> result3 = implicit_euler(funcs, 0.0, u0, tend, 0.05);
	exportData(result, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test1/implicit_euler.txt");
	for (const auto & i : result3)
	{
		std::cout << i << std::endl;
	}
	std::vector<std::vector<double>> result4 = trapezoidal_rule_method(funcs, 0.0, u0, tend, 0.05);
	exportData(result, "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab1 (Numerical methods for solving ordinary differential equations)/data/test1/symmetric.txt");
	for (const auto & i : result4)
	{
		std::cout << i << std::endl;
	}


}