#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <memory>
#include "../../Lab1/src/matrix_operations.h"
#include "../../Lab5/src/FindRoots.h"


int main()
{
	auto ex1 = [](double x) { return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75); };
	auto dex1 = [](double x) { return 0.643225 - 4.6844 * x + 11.9475 * x * x - 12.88 * x * x * x + 5 * x * x * x * x; };
	double a1=0; double b1=1;
	auto ex2 = [](double x) {return std::sqrt(x + 1) - 1;};
	auto dex2 = [](double x) {return 1.0 / (2 * std::sqrt(1 + x));};
	double a2=-1; double b2=10;
	auto ex3 = [](double x) {return 35 * x * x * x - 67 * x * x - 3 * x + 3;};
	auto dex3 = [](double x) {return -3 - 134 * x + 105 * x * x;};
	double a3=0; double b3=1;
	int initialGridPoints = 15;
	double maxPoints = 2000;
	double tolerance = 10e-6;
	std::cout << "[bisectionMethod] ex1: " << bisectionMethod<double>(ex1, a1, b1, initialGridPoints, tolerance) << '\n';
	std::cout << "[bisectionMethod] ex2: " << bisectionMethod<double>(ex2, a2, b2, initialGridPoints, tolerance) << '\n';
	std::cout << "[bisectionMethod] ex3: " << bisectionMethod<double>(ex3, a3, b3, initialGridPoints, tolerance) << '\n';
	std::cout << '\n';
	std::cout << "[newtonMethod (analytical derivative)] ex1: " << newtonMethod<double>(ex1, a1, b1, tolerance, 1000, dex1) << '\n';
	std::cout << "[newtonMethod (analytical derivative)] ex2: " << newtonMethod<double>(ex2, a2, b2, tolerance, 1000, dex2) << '\n';
	std::cout << "[newtonMethod (analytical derivative)] ex3: " << newtonMethod<double>(ex3, a3, b3, tolerance, 1000, dex3) << '\n';
	std::cout << '\n';
	std::cout << "[newtonMethod (numerical derivative)] ex1: " << newtonMethod<double>(ex1, a1, b1, tolerance, 1000) << '\n';
	std::cout << "[newtonMethod (numerical derivative)] ex2: " << newtonMethod<double>(ex2, a2, b2, tolerance, 1000) << '\n';
	std::cout << "[newtonMethod (numerical derivative)] ex3: " << newtonMethod<double>(ex3, a3, b3, tolerance, 1000) << '\n';

}
