#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <fstream>

#include "../../Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../Lab3 (Solving interpolation problems)/src/grid_generators.h"

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
		T x = 0;
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

template <typename T>
std::pair<T, int> newtonMethod(std::function<T(T)> func, std::function<T(T)> derivative, T x0, T tolerance, int maxIterations) {
	T x = x0;
	int iterations = 0;

	for (int i = 0; i < maxIterations; ++i) {
		T f = func(x);
		T df = derivative(x);

//		if (std::abs(df) < tolerance) {
//			std::cerr << "Производная близка к нулю. Метод может не сойтись.\n";
//			break;
//		}

		T xNext = x - f / df;

		if (std::abs(xNext - x) < tolerance) {
			return {xNext, iterations + 1};
		}

		x = xNext;
		iterations++;
	}

	return {x, iterations};
}

template <typename T>
bool isPointAlreadyFound(const std::vector<std::vector<T>>& results, const std::vector<T>& point, T tolerance) {
	for (const auto& result : results) {
		if (compute2Norm(result - point) < tolerance) {
			return true;
		}
	}
	return false;
}

template <typename T>
bool isPointAlreadyFound(const std::vector<std::pair<std::vector<T>, int>>& results, const std::vector<T>& point, T tolerance) {
	for (const auto& result : results) {
		if (compute2Norm(result.first - point) < tolerance) {
			return true;
		}
	}
	return false;
}

template <typename T>
std::vector<std::pair<std::vector<T>, int>> newtonMethod(
	const std::function<std::vector<T>(const std::vector<T>&)>& funcs,
	T tolerance, int maxIterations,
	const std::function<std::vector<T>(const std::pair<std::vector<std::vector<T>>, std::vector<T>>&)>& solve,
	T ax, T bx, T nx, T ay, T by, T ny
) {


	T h = tolerance;
	auto grid = generate_uniform_grid(ax, bx, nx, ay, by, ny);
	std::vector<std::pair<std::vector<T>, int>> results;
	results.reserve(grid.size());

	std::ofstream outputFileStream("/Users/ivandybko/Projects/Numerical_methods/Lab5/src/data/sys2_1.txt");
	if (!outputFileStream.is_open()) {
		std::cerr << "Ошибка при открытии файла для записи!" << std::endl;
		return results;
	}
	auto compute_jacobian = [&](const std::vector<T>& x) {
	  size_t n = x.size();
	  std::vector<std::vector<T>> J(n, std::vector<T>(n));
	  for (size_t i = 0; i < n; ++i) {
		  std::vector<T> x_plus = x, x_minus = x;
		  for (size_t j = 0; j < n; ++j) {
			  x_plus[j] += h;
			  x_minus[j] -= h;
			  J[i][j] = (funcs(x_plus)[i] - funcs(x_minus)[i]) / (2 * h);
			  x_plus[j] = x[j];
			  x_minus[j] = x[j];
		  }
	  }
	  return J;
	};
	for (const auto& point : grid) {
		std::vector<T> x = {point.first, point.second};
		int iterations = 0;

		for (; iterations < maxIterations; ++iterations) {
			std::vector<T> f = funcs(x);
			if (compute2Norm(f) < tolerance) {
				break;
			}
			std::vector<std::vector<T>> J = compute_jacobian(x);
			std::vector<T> dx = solve({J, f});
			if (dx.empty()){break;};
			for (auto& dxi : dx) {
				dxi = -dxi;
			}
			x=x+dx;
			if (!(x[0] >= ax && x[0] <= bx) or !(x[1] >= ay && x[1] <= by)){
				iterations=maxIterations+10;
				break;

			}
			if (compute2Norm(dx) < tolerance) {
				if (!isPointAlreadyFound(results,x,tolerance)){
					results.emplace_back(x, iterations);
				}
				break;
			}
		}
		int finalIterations = (iterations >= maxIterations) ? -1 : iterations;

		outputFileStream << point.first << " " << point.second << " " << finalIterations << "\n";
	}
	outputFileStream.close();
	return results;
}

template <typename T>
std::vector<std::vector<T>> newtonMethod(
    const std::function<std::vector<T>(const std::vector<T>&)>& funcs,
    T tolerance,
    int maxIterations,
    const std::function<std::vector<T>(const std::pair<std::vector<std::vector<T>>, std::vector<T>>&)>& solve,
    const std::vector<std::pair<T, T>>& bounds,
    const std::vector<size_t>& subdivisions
) {
    T h = tolerance;
    auto grid = generate_nd_uniform_grid(bounds, subdivisions);

    std::vector<std::vector<T>> results;
    results.reserve(grid.size());

    auto computeJacobian = [&](const std::vector<T>& x) -> std::vector<std::vector<T>> {
        size_t n = x.size();
        std::vector<std::vector<T>> J(n, std::vector<T>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                auto x_plus  = x;
                auto x_minus = x;
                x_plus[j]  += h;
                x_minus[j] -= h;
                J[i][j] = (funcs(x_plus)[i] - funcs(x_minus)[i]) / (2 * h);
            }
        }
        return J;
    };

    for (const auto& point : grid) {
        std::vector<T> x = point;
        int iterations = 0;
        for (; iterations < maxIterations; ++iterations) {
            auto f_val = funcs(x);
            if (compute2Norm(f_val) < tolerance) {
                break;
            }
            auto J = computeJacobian(x);
            auto dx = solve({ J, f_val });
            if (dx.empty()) {
                break;
            }
            for (auto& dxi : dx) {
                dxi = -dxi;
            }
            x = x + dx;
            bool inBounds = true;
            for (size_t i = 0; i < x.size(); ++i) {
                if (x[i] < bounds[i].first || x[i] > bounds[i].second) {
                    inBounds = false;
                    break;
                }
            }
            if (!inBounds) {
                iterations = maxIterations + 10;
                break;
            }
        	f_val = funcs(x);
        	if (compute2Norm(f_val) < tolerance || compute2Norm(dx) < tolerance) {
        		if (!isPointAlreadyFound(results, x, tolerance)) {
        			results.emplace_back(x);
        		}
        		break;
        	}
        }
    }
    return results;
}


