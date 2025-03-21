#ifndef LAB3_LAB1_SRC_DATA_TRIDIAGONALMATRIXALGORITHM_H_
#define LAB3_LAB1_SRC_DATA_TRIDIAGONALMATRIXALGORITHM_H_
#include <vector>
#include <stdexcept>

template <typename T>
std::vector<T> TridiagonalMatrixAlgorithm(const std::vector<T>& a,
	const std::vector<T>& b,
	const std::vector<T>& c,
	const std::vector<T>& d) {
	size_t n = b.size();
	if (a.size() != n - 1 || c.size() != n - 1 || d.size() != n) {
		throw std::invalid_argument("Размеры диагоналей не соответствуют размеру системы.");
	}
	std::vector<T> cp(n - 1);
	std::vector<T> dp(n);
	cp[0] = c[0] / b[0];
	dp[0] = d[0] / b[0];
	for (size_t i = 1; i < n; ++i) {
		T denom = b[i] - a[i - 1] * cp[i - 1];
		if (denom == 0.0) {
			throw std::runtime_error("Система вырождена или плохо обусловлена.");
		}
		if (i < n - 1) {
			cp[i] = c[i] / denom;
		}
		dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / denom;
	}
	std::vector<T> x(n);
	x[n - 1] = dp[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		x[i] = dp[i] - cp[i] * x[i + 1];
	}
	return x;
}
#endif //LAB3_LAB1_SRC_DATA_TRIDIAGONALMATRIXALGORITHM_H_
