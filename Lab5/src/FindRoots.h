#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#include "../../Lab3/src/grid_generators.h"

template <typename T>
std::vector<std::pair<T, T>> rootIntervals(
	std::function<T(T)> func, T a, T b, int initialGridPoints, T maxPoints) {
	std::vector<std::pair<T, T>> intervals;
	auto grid = generate_uniform_grid(a, b, initialGridPoints, func);

	for (size_t i = 1; i < grid->size(); ++i) {
		T x1 = (*grid)[i - 1].first;
		T f1 = (*grid)[i - 1].second;
		T x2 = (*grid)[i].first;
		T f2 = (*grid)[i].second;
		if (f1 * f2 <= 0) {
			intervals.emplace_back(x1, x2);
		}
	}
	if (intervals.empty() && initialGridPoints < maxPoints) {
		return rootIntervals(func, a, b, initialGridPoints * 2, maxPoints);
	}
	return intervals;
}

template <typename T>
std::vector<std::pair<T, int>> bisectionMethod(
	std::function<T(T)> func, T a, T b, int initialGridPoints, T tolerance) {
	auto intervals = rootIntervals(func, a, b, initialGridPoints, tolerance);
	std::vector<std::pair<T, int>> roots;
	for (const auto& interval : intervals) {
		T left = interval.first;
		T right = interval.second;
		int iterations = 0;
		while (right - left > tolerance) {
			T mid = (left + right) / 2;
			if (func(left) * func(mid) <= 0) {
				right = mid;
			} else {
				left = mid;
			}
			iterations++;
		}
		roots.emplace_back((left + right) / 2, iterations);
	}
	return roots;
}

template <typename T>
std::vector<std::pair<T, int>> newtonMethod(std::function<T(T)> func, T a, T b, T tolerance, int maxIterations, std::function<T(T)> derivative = nullptr) {
	int initialGridPoints = (b-a)*20; T h = tolerance;
	auto intervals = rootIntervals(func, a, b, initialGridPoints, tolerance);
	std::vector<std::pair<T, int>> roots;
	for (const auto& interval : intervals) {
		int iterationsForCurrentRoot = 0;
		T x0 = (interval.first + interval.second) / 2;
		T x = x0;
		for (int i = 0; i < maxIterations; ++i) {
			T f = func(x);
			T df;
			if (derivative) {
				df = derivative(x);
			} else {
				df = (func(x + h) - func(x)) / (h);
			}
			if (std::abs(df) < tolerance) {
				std::cerr << "Производная близка к нулю на интервале [" << interval.first << ", " << interval.second << "].\n";
				break;
			}
			T xNext = x - f / df;

			if (xNext < interval.first) {
				xNext = interval.first;
			} else if (xNext > interval.second) {
				xNext = interval.second;
			}

			if (std::abs(xNext - x) < tolerance) {
				roots.push_back({xNext, iterationsForCurrentRoot});
				break;
			}
			x = xNext;
			iterationsForCurrentRoot++;
			if (i == maxIterations - 1) {
				std::cerr << "Превышено максимальное количество итераций для интервала [" << interval.first << ", " << interval.second << "].\n";
			}
		}
	}
	return roots;
}

