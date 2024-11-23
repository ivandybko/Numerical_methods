#include <iostream>
#include <fstream>
#include "grid_generators.h"
#include "NewtonInterpolation.h"
#include "../../Lab1/src/matrix_operations.h"
template <typename T>
void exportGridToFile(const std::unique_ptr<std::vector<std::pair<T, T>>>& grid, const std::string& filename) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Ошибка открытия файла для записи!" << std::endl;
		return;
	}
	for (const auto& point : *grid) {
		outFile << point.first << " " << point.second << std::endl;
	}
	outFile.close();

}

template <typename T>
T calculateErrorNorms(
	const std::function<T(T)>& P,
	const std::vector<std::pair<T, T>>& testGrid) {
	T maxError = 0;

	for (auto xi : testGrid) {
		T error = std::abs(xi.second - P(xi.first));
		maxError = std::max(maxError, error);
	}
	return maxError;
}


int main()
{
	double a = -1.0, b = 1.0;
	std::function<double(double)> func = [](double x) { return 1/(1+25*x*x); };
	for (size_t n: {4,8,16,32})
	{
		auto uniformPoints = generate_uniform_grid(a, b, n, func);
		auto uniformPointsTest = generate_uniform_grid(a, b, 1000, func);
		auto uniformPolynomial = NewtonInterpolation(*uniformPoints);
		exportGridToFile(generate_uniform_grid(a, b, 1000, uniformPolynomial),"/Users/ivandybko/Projects/Numerical_methods/Lab3/src/data/Runge/Runge_"+std::to_string(n)+"_uniform_grid.txt");
		std::cout << "Норма ошибки интерполяции Лагранжа на равномерной сетке при n=" << n  << ": " << calculateErrorNorms<double>(uniformPolynomial, *uniformPointsTest) << '\n';
		auto chebyshevPoints = generate_chebyshev_grid(a, b, n, func);
		auto chebyshevPointsTest = generate_chebyshev_grid(a, b, 1000, func);
		auto chebyshevPolynomial = NewtonInterpolation(*chebyshevPoints);
		exportGridToFile(generate_chebyshev_grid(a, b, 1000, chebyshevPolynomial),"/Users/ivandybko/Projects/Numerical_methods/Lab3/src/data/Runge/Runge_"+std::to_string(n)+"_chebyshev_grid.txt");
		std::cout << "Норма ошибки интерполяции Лагранжа на Чебышевской сетке при n=" << n  << ": " << calculateErrorNorms<double>(chebyshevPolynomial, *chebyshevPointsTest) << "\n\n" ;
	}
}
