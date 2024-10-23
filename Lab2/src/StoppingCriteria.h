#ifndef LAB_2_SRC_STOPPINGCRITERIA_H_
#define LAB_2_SRC_STOPPINGCRITERIA_H_
#include <iostream>
#include <vector>

template <typename T>
bool checkConvergence(const std::vector<T> &x, const std::vector<T> &x_old, T eps, T norm_C = 0.5)
{
	size_t n = x.size();
	T criteria{ 0.0 };
	for (size_t i = 0; i < n; ++i)
	{
		criteria += abs(x[i] - x_old[i]);
	}
	return criteria < (1-norm_C)/(norm_C) * eps;
}


template <typename T>
bool checkConvergence(const std::vector<T> & Ax_minus_f, T eps)
{
	size_t n = Ax_minus_f.size();
	T criteria{ 0.0 }, criteria2{ 0.0 };
	for (size_t i = 0; i < n; ++i)
	{
		criteria += abs(Ax_minus_f[i]);
	}
	return criteria < eps;
}
#endif //LAB_2_SRC_STOPPINGCRITERIA_H_
