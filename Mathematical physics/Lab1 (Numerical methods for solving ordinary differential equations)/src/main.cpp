#include <iostream>
#include <vector>
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "runge_kutta_4th_order.h"
int main()
{
	std::vector<std::function<double(double, const std::vector<double>&)>> funcs = {
		[](double t, const std::vector<double>& u) { return u[1]  ; },
		[](double t, const std::vector<double>& u) { return  - u[0] ; }
	};

	std::vector<double> u0 = {1.0, 0.0};
	std::vector<std::vector<double>> result = RungeKutta4System(funcs, 0.0, u0, 10, 0.1);
	for (const auto & i : result)
	{
		std::cout << i << std::endl;
	}
}