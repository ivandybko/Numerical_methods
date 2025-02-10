#ifndef LAB3_SRC_NEWTONINTERPOLATION_H_
#define LAB3_SRC_NEWTONINTERPOLATION_H_
#include <iostream>
#include <vector>
#include <functional>

template <typename T>
std::function<T(T)> NewtonInterpolation(const std::vector<std::pair<T, T>>& points){
	if (points.empty()) {
		throw std::invalid_argument("Вектор точек пуст.");
	}
	size_t n = points.size();
	std::vector<std::vector<T>> dividedDifferences(n, std::vector<T>(n, 0));
	for (size_t i = 0; i < n; ++i) {
		dividedDifferences[i][0] = points[i].second;
	}
	for (size_t j = 1; j < n; ++j) {
		for (size_t i = 0; i < n - j; ++i) {
			dividedDifferences[i][j] =
				(dividedDifferences[i + 1][j - 1] - dividedDifferences[i][j - 1]) /
					(points[i + j].first - points[i].first);
		}
	}

	return [points, dividedDifferences, n](T x) {
	  T result = dividedDifferences[0][0];
	  T product = 1;

	  for (size_t i = 1; i < n; ++i) {
		  product *= (x - points[i - 1].first);
		  result += product * dividedDifferences[0][i];
	  }

	  return result;
	};
};
#endif //LAB3_SRC_NEWTONINTERPOLATION_H_
