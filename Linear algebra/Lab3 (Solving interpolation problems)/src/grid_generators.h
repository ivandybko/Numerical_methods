#ifndef LAB3_SRC_GRID_GENERATORS_H_
#define LAB3_SRC_GRID_GENERATORS_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <functional>

template <typename T>
std::unique_ptr<std::vector<std::pair<T, T>>> generate_uniform_grid(T a, T b, int n, std::function<T(T)> func){
	if (n < 2) {
		throw std::invalid_argument("Количество точек должно быть >= 2.");
	}
	if (a >= b) {
		throw std::invalid_argument("Левая граница должна быть меньше правой.");
	}
	auto points = std::make_unique<std::vector<std::pair<T, T>>>();
	T h = (b - a) / (n - 1);
	for (size_t i = 0; i < n; ++i) {
		T x = a + i * h;
		points->push_back({x, func(x)});
	}
	return points;
};

template <typename T>
std::unique_ptr<std::vector<std::pair<T, T>>> generate_chebyshev_grid(T a, T b, int n, std::function<T(T)> func){
	if (n < 2) {
		throw std::invalid_argument("Количество точек должно быть >= 2.");
	}
	if (a >= b) {
		throw std::invalid_argument("Левая граница должна быть меньше правой.");
	}
	auto points = std::make_unique<std::vector<std::pair<T, T>>>();
	for (size_t i = 0; i < n; ++i) {
		T x = 0.5 * (a + b) + 0.5 * (b - a) * cos(M_PI * (2 * i + 1) / (2 * (n+1)));
		points->push_back({x, func(x)});
	}
	std::reverse(points->begin(), points->end());
	return points;
};

template <typename T>
std::vector<std::pair<T, T>> generate_uniform_grid(
	T ax, T bx, int nx, T ay, T by, int ny) {
	if (nx < 2 || ny < 2) {
		throw std::invalid_argument("Количество точек должно быть >= 2 в обоих направлениях.");
	}
	if (ax >= bx || ay >= by) {
		throw std::invalid_argument("Левая граница должна быть меньше правой в обоих направлениях.");
	}
	std::vector<std::pair<T, T>> points;
	points.reserve(nx * ny);

	T hx = (bx - ax) / (nx);
	T hy = (by - ay) / (ny);

	for (int i = 0; i < nx; ++i) {
		T x = ax + i * hx;
		for (int j = 0; j < ny; ++j) {
			T y = ay + j * hy;
			points.emplace_back(x, y);
		}
	}

	return points;
}

template <typename T>
std::vector<std::vector<T>> generate_nd_uniform_grid(
	const std::vector<std::pair<T, T>>& bounds,
	const std::vector<size_t>& subdivisions
) {
	size_t dims = bounds.size();
	size_t totalPoints = 1;
	std::vector<T> steps(dims);

	for (size_t d = 0; d < dims; ++d) {
		totalPoints *= subdivisions[d];
		if (subdivisions[d] > 1)
			steps[d] = (bounds[d].second - bounds[d].first) / static_cast<T>(subdivisions[d] - 1);
		else
			steps[d] = T{};
	}

	std::vector<size_t> strides(dims, 1);
	for (size_t d = 1; d < dims; ++d) {
		strides[d] = strides[d - 1] * subdivisions[d - 1];
	}

	std::vector<std::vector<T>> grid;
	grid.reserve(totalPoints);

	for (size_t index = 0; index < totalPoints; ++index) {
		std::vector<T> point(dims);
		for (size_t d = 0; d < dims; ++d) {
			size_t coordIndex = (index / strides[d]) % subdivisions[d];
			point[d] = bounds[d].first + static_cast<T>(coordIndex) * steps[d];
		}
		grid.push_back(std::move(point));
	}

	return std::move(grid);
}
#endif //LAB3_SRC_GRID_GENERATORS_H_
