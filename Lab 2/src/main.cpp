#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include "SimpleIterationMethod.h"
#include "JacobiMethod.h"
#include "SeidelMethod.h"
#include "SuccessiveRelaxationMethod.h"

template <typename T>
std::unique_ptr<std::pair<std::vector<std::vector<T>>, std::vector<T>>> readData(const std::string& path) {
	auto data = std::make_unique<std::pair<std::vector<std::vector<T>>, std::vector<T>>>();
	std::vector<std::vector<T>>& matrix = data->first;
	std::vector<T>& vector = data->second;
	std::ifstream inputMatrix(path + "/matrix.dat");
	if (!inputMatrix.is_open()) {
		throw std::invalid_argument("Файл matrix.dat не найден");
		exit(1);
	}
	int n;
	inputMatrix >> n;
	matrix.resize(n);
	for (auto& row : matrix) {
		row.resize(n);
	}
	for (auto& row : matrix) {
		for (auto& element : row) {
			inputMatrix >> element;
		}
	}
	inputMatrix.close();
	std::ifstream inputVector(path + "/vector.dat");
	if (!inputVector.is_open()) {
		throw std::invalid_argument("Файл vector.dat не найден");
		exit(1);
	}
	inputVector >> n;
	vector.resize(n);
	for (auto& element : vector) {
		inputVector >> element;
	}
	inputVector.close();
	return data;
}

template <typename T>
std::tuple<std::shared_ptr<std::vector<T>>, std::shared_ptr<std::vector<T>>,
		   std::shared_ptr<std::vector<T>>, std::shared_ptr<std::vector<T>>>
TridiagonalMatrixGenerator(int N) {
	int n = 200 + N;
	auto a = std::make_shared<std::vector<T>>(n-1, 1);
	auto b = std::make_shared<std::vector<T>>(n, 4);
	auto c = std::make_shared<std::vector<T>>(n-1, 1);
	auto d = std::make_shared<std::vector<T>>(n);
	(*d)[0] = 6;
	for (int i = 1; i < n - 1; ++i) {
		(*d)[i] = 10 - 2 * ((i+1) % 2);
	}
	(*d)[n - 1] = 9 - 3 * (n % 2);
	return std::make_tuple(a, b, c, d);
}

int main()
{
	std::cout << "Введите путь к директории: ";
	std::string path;
	std::cin >> path;
	std::cout << '\n';
	auto data = *readData<double>(path);
	auto sol_simpleiter = SimpleIterationMethod<double>(data, 100000,1e-2,1e-7);
	if (!sol_simpleiter.x.empty()){std::cout << "Решение найденное за " << sol_simpleiter.iter << " итераций методом простой итерации (норма матрицы С=" << sol_simpleiter.norm_C << "): " << sol_simpleiter.x << '\n';}
	auto sol_jacobi = JacobiMethod<double>(data, 100000,1e-2);
	if (!sol_jacobi.x.empty()){std::cout << "Решение найденное за " << sol_jacobi.iter << " итераций методом Якоби (норма матрицы С=" << sol_jacobi.norm_C << "): " << sol_jacobi.x << '\n';}
	auto sol_seidel = SeidelMethod<double>(data, 100000,1e-2);
	if (!sol_seidel.x.empty()){std::cout << "Решение найденное за " << sol_seidel.iter << " итераций методом Зейделя (сумма норм матриц ||С_L||+||C_U||=" << sol_seidel.norm_C << "): " << sol_seidel.x << '\n';}
	auto sol_relaxation = SuccessiveRelaxationMethod<double>(data, 10000,1e-2,0.4);
	if (!sol_relaxation.x.empty()){std::cout << "Решение найденное за " << sol_relaxation.iter << " итераций методом релаксации (сумма норм матриц ||С_L||+||C_U||=" << sol_relaxation.norm_C << "): " << sol_relaxation.x << '\n';}
	auto tridiagonal_matrix = TridiagonalMatrixGenerator<double>(7);
	auto sol_seidel_tridiag = SeidelMethod<double>(*std::get<0>(tridiagonal_matrix),
															  *std::get<1>(tridiagonal_matrix),
			  												  *std::get<2>(tridiagonal_matrix),
			  												  *std::get<3>(tridiagonal_matrix), 10000,1e-4);
	if (!sol_seidel_tridiag.x.empty()){std::cout << "Решение найденное за " << sol_seidel_tridiag.iter << " итераций методом Зейделя случая трехдиагональной матрицы: "<< sol_seidel_tridiag.x << '\n';}
	auto sol_relaxation_tridiag = SuccessiveRelaxationMethod<double>(*std::get<0>(tridiagonal_matrix),
		*std::get<1>(tridiagonal_matrix),
		*std::get<2>(tridiagonal_matrix),
		*std::get<3>(tridiagonal_matrix), 10000,1e-4,1.1);
	if (!sol_relaxation_tridiag.x.empty()){std::cout << "Решение найденное за " << sol_relaxation_tridiag.iter << " итераций методом Зейделя случая трехдиагональной матрицы: "<< sol_relaxation_tridiag.x << '\n';}
}
