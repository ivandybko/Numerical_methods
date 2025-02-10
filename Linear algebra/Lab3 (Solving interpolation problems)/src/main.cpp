#include <iostream>
#include <fstream>
#include <memory>
#include "grid_generators.h"
#include "NewtonInterpolation.h"
#include "../../Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "SplineInterpolation.h"
template <typename T>
void exportGridToFile(const std::unique_ptr<std::vector<std::pair<T, T>>>& grid, const std::string& filename) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Ошибка открытия файла для записи!" << std::endl;
		return;
	}
	for (const auto& point : *grid) {
		outFile << std::setprecision(40) << point.first << " " << std::setprecision(18) << point.second << std::endl;
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
	std::function<double(double)> func = [](double x) { return x*x+x; };
	for (size_t n : {64,128,256,512,1024}){
		auto uniformPoints = generate_uniform_grid(a, b, n, func);
		CubicSpline<double> spline_x;
		spline_x.build(*uniformPoints,2,2);
		std::function<double(double)> interpolator = [&spline_x](double x)
		{
		  return spline_x.interpolate(x);
		};
		auto uniformPointsTest = generate_uniform_grid(a, b, 10000, func);
		std::cout << "Норма ошибки сплайн-интерполяции при n=" << n << ": "
				  << calculateErrorNorms<double>(interpolator, *uniformPointsTest) << '\n';
//		auto z1 = calculateErrorNorms<double>(interpolator, *uniformPointsTest)/err_old_u;
//		std::cout << z1 << ' '  << "log" << 0.5 << "("") = " << std::log(z1) / std::log(0.5) << std::endl;
//		err_old_u = calculateErrorNorms<double>(interpolator, *uniformPointsTest);
//		exportGridToFile(generate_uniform_grid(a, b, 10000000, interpolator),
//			"/Users/ivandybko/Projects/Numerical_methods/Lab3/src/data/x^2/x^2_"+std::to_string(n)+"_spline_uniform_grid.txt");
//		auto chebyshevPoints = generate_chebyshev_grid(a, b, n, func);
//		CubicSpline<double> spline_x_cheb;
//		spline_x_cheb.build(*chebyshevPoints);
//		std::function<double(double)> interpolator_cheb = [&spline_x_cheb](double x)
//		{
//		  return spline_x_cheb.interpolate(x);
//		};
//		auto a_real = (*chebyshevPoints)[0].first;
//		auto b_real = (*chebyshevPoints)[n-1].first;
//		auto chebyshevPointsTest = generate_chebyshev_grid(a_real, b_real, 100, func);
//		std::cout << "Норма ошибки сплайн-интерполяции на Чебышевской сетке при n=" << n << ": "
//				  << calculateErrorNorms<double>(interpolator_cheb, *chebyshevPointsTest) << '\n';
//		exportGridToFile(generate_chebyshev_grid(a_real, b_real, 1000, interpolator_cheb),"/Users/ivandybko/Projects/Numerical_methods/Lab3/src/data/x^2/x^2_"+std::to_string(n)+"_spline_chebyshev_grid.txt");
	}
}
