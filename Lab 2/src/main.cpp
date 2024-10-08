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

int main()
{
	std::cout << "Введите путь к директории: ";
	std::string path;
	std::cin >> path;
	std::cout << '\n';
	auto data = *readData<double>(path);
	auto x_simpleiter = SimpleIterationMethod<double>(data, 10000,1e-8);
	std::cout << "Решение найденное методом простой итерации: " << x_simpleiter << '\n';
	auto x_jacobi = JacobiMethod<double>(data, 10000,1e-2);
	std::cout << "Решение найденное методом Якоби: "<< x_jacobi << '\n';
	auto x_seidel = SeidelMethod<double>(data, 10000,1e-2);
	std::cout << "Решение найденное методом Зейделя: "<< x_seidel << '\n';
	auto x_relaxation = SuccessiveRelaxationMethod<double>(data, 10000,1e-2,1.2);
	std::cout << "Решение найденное методом релаксации: "<< x_relaxation << '\n';
}
