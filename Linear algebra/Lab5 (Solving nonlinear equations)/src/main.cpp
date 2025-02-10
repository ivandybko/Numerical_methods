#include <iostream>
#include <cmath>
#include "../../Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../Lab5 (Solving nonlinear equations)/src/FindRoots.h"
#include "../../Lab1 (Direct methods for solving systems of linear equations)/src/Gauss.h"


//int main()
//{
//	std::cout << "Singular nonlinear equations: \n";
//	auto ex1 = [](double x) { return (x - 1) * (x - 1) * (x - 1) ; };
//	auto dex1 = [](double x) { return 0.643225 - 4.6844 * x + 11.9475 * x * x - 12.88 * x * x * x + 5 * x * x * x * x; };
//	double a1=0, b1=1;
//	auto ex2 = [](double x) {return std::sqrt(x + 1) - 1;};
//	auto dex2 = [](double x) {return 1.0 / (2 * std::sqrt(1 + x));};
//	double a2=-1, b2=10;
//	auto ex3 = [](double x) {return 35 * x * x * x - 67 * x * x - 3 * x + 3;};
//	auto dex3 = [](double x) {return -3 - 134 * x + 105 * x * x;};
//	double a3=0, b3=1;
//	int initialGridPoints = 15;
//	double maxPoints = 2000;
//	double tolerance = 10e-6;
////	std::cout << "[bisectionMethod] ex1: " << bisectionMethod<double>(ex1, a1, b1, initialGridPoints, tolerance) << '\n';
////	std::cout << "[bisectionMethod] ex2: " << bisectionMethod<double>(ex2, a2, b2, initialGridPoints, tolerance) << '\n';
////	std::cout << "[bisectionMethod] ex3: " << bisectionMethod<double>(ex3, a3, b3, initialGridPoints, tolerance) << '\n';
////	std::cout << '\n';
////	std::cout << "[newtonMethod (analytical derivative)] ex1: " << newtonMethod<double>(ex1, a1, b1, tolerance, 1000, dex1) << '\n';
////	std::cout << "[newtonMethod (analytical derivative)] ex2: " << newtonMethod<double>(ex2, a2, b2, tolerance, 1000, dex2) << '\n';
////	std::cout << "[newtonMethod (analytical derivative)] ex3: " << newtonMethod<double>(ex3, a3, b3, tolerance, 1000, dex3) << '\n';
////	std::cout << '\n';
//	std::cout << "[newtonMethod (numerical derivative)] ex1: " << newtonMethod<double>(ex1, a1, b1, tolerance, 1000) << '\n';
////	std::cout << "[newtonMethod (numerical derivative)] ex2: " << newtonMethod<double>(ex2, a2, b2, tolerance, 1000) << '\n';
////	std::cout << "[newtonMethod (numerical derivative)] ex3: " << newtonMethod<double>(ex3, a3, b3, tolerance, 1000) << '\n';
//	std::cout << "Systems of nonlinear equations: \n";
//	auto sys1 = [](const std::vector<double>& x) {
//	  return std::vector<double>{
//		  x[0] * x[0] - x[1] * x[1] - 15,
//		  x[0] * x[1] + 4
//	  };
//	};
////	double ax1=-10,bx1=10, ay1=-10,by1=10;
////	std::cout << "[newtonMethod] sys1: " << newtonMethod<double>(sys1,tolerance, 30, Gauss<double>, ax1, bx1, (std::abs(ax1)+std::abs(bx1))*10, ay1, by1, (std::abs(ay1)+std::abs(by1))*10) << '\n';
////	auto sys2 = [](const std::vector<double>& x) {
////	  return std::vector<double>{
////		  x[0] * x[0] + x[1]*x[1] + x[0] + x[1] - 8,
////		  x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 7
////	  };
////	};
////	std::cout << "[newtonMethod] sys2: " << newtonMethod<double>(sys2,tolerance, 30, Gauss<double>, ax1, bx1, abs(ax1*bx1)*2, ay1, by1, abs(ay1*by1)*2) << '\n';
////
//
//}

int main() {
	auto func = [](double x) { return std::pow(x - 1, 3); };
	auto derivative = [](double x) { return 3 * std::pow(x - 1, 2); };

	double x0 = 0.0;

	double tolerance = 1e-6;
	int maxIterations = 100;

	auto [root, iterations] = newtonMethod<double>(func, derivative, x0, tolerance, maxIterations);

	std::cout << "Найденный корень: " << root << "\n";
	std::cout << "Число итераций: " << iterations << "\n";

	return 0;
}
