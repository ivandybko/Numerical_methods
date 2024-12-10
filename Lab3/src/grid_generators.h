#ifndef LAB3_SRC_GRID_GENERATORS_H_
#define LAB3_SRC_GRID_GENERATORS_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>

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

#endif //LAB3_SRC_GRID_GENERATORS_H_
