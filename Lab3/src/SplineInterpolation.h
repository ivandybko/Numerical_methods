#ifndef LAB3_SRC_SPLINEINTERPOLATION_H_
#define LAB3_SRC_SPLINEINTERPOLATION_H_
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include "../../Lab1/src/TridiagonalMatrixAlgorithm.h"

template <typename T>
class CubicSpline {
 private:
	std::vector<T> a, b, c, d; // Коэффициенты сплайна
	std::vector<T> x;         // Узлы сетки

 public:
	void build(const std::vector<std::pair<T, T>>& grid) {
		size_t n = grid.size();
		if (n < 2) {
			throw std::invalid_argument("Grid must contain at least two points.");
		}
		x.resize(n);
		a.resize(n);
		for (size_t i = 0; i < n; ++i) {
			x[i] = grid[i].first;
			a[i] = grid[i].second;
		}
		std::vector<T> h(n - 1);
		for (size_t i = 0; i < n - 1; ++i) {
			h[i] = x[i + 1] - x[i];
			if (h[i] <= 0) {
				throw std::invalid_argument("Grid points must be sorted and unique.");
			}
		}
		std::vector<T> sub_diag(n - 2);
		std::vector<T> main_diag(n - 1);
		std::vector<T> sup_diag(n - 2);
		std::vector<T> rhs(n - 1);
		main_diag[0] = 1;
		main_diag[n - 1] = 1;
		rhs[0] = 0;
		rhs[n - 1] = 0;
		for (size_t i = 1; i < n - 1; ++i) {
			sub_diag[i - 1] = h[i - 1];
			main_diag[i] = 2 * (h[i - 1] + h[i]);
			sup_diag[i - 1] = h[i];
			rhs[i] = 3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]);
		}
		c.resize(n);
		c = TridiagonalMatrixAlgorithm<T>(sub_diag, main_diag, sup_diag, rhs);
		c[0] = c[n - 1] = 0;
		b.resize(n - 1);
		d.resize(n - 1);
		for (size_t i = 0; i < n - 1; ++i) {
			b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
			d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
		}
	}

	T interpolate(T point) const {
		if (x.empty()) {
			throw std::runtime_error("Spline has not been built.");
		}
		size_t i = 0;
		for (; i < x.size() - 1; ++i) {
			if (point <= x[i + 1]) {
				break;
			}
		}
		if (i >= x.size() - 1) {
			throw std::out_of_range("Point is out of the interpolation range.");
		}
		T dx = point - x[i];
		return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
	}
};
#endif //LAB3_SRC_SPLINEINTERPOLATION_H_