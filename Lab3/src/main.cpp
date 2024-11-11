#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "FrancisKublanovskayaQR.h"
#include "InverseIterationMethod.h"
template <typename T>
std::unique_ptr<std::vector<std::vector<T>>> readData(const std::string& path) {
	auto data = std::make_unique<std::vector<std::vector<T>>>();
	std::ifstream inputMatrix(path);
	if (!inputMatrix.is_open()) {
		throw std::invalid_argument("Файл не найден");
	}
	int n;
	inputMatrix >> n;
	data->resize(n);
	for (auto& row : *data) {
		row.resize(n);
	}
	for (auto& row : *data) {
		for (auto& element : row) {
			inputMatrix >> element;
		}
	}
	inputMatrix.close();
	return data;
}


int main()
{
	std::cout << "Введите путь к директории: ";
	std::string path;
	std::cin >> path;
	std::cout << '\n';
	auto data = readData<double>(path);
	auto eigen_values = QRAlgorithm<double>(*data, 1e-10, true, true);
	std::cout << eigen_values << '\n';
	std::cout << "Решение методом обратной итерации без соотношения Рэлея: \n";
	auto eigen_vectors = InverseIterationMethod<double>(*data, eigen_values, 1e-10);
	for (int i = 0; i < eigen_values.size(); ++i)
	{
		std::cout << "λ=" << eigen_values[i] << "; Eigen vector: " << *eigen_vectors[i] <<";\n";
	}
	std::cout << "Решение методом обратной итерации с соотношением Рэлея: \n";
	auto eigen_system = InverseIterationMethod<double>(*data, 1e-10);
	for (int i = 0; i < eigen_system.first.size(); ++i)
	{
		std::cout << "λ=" << eigen_system.first[i] << "; Eigen vector: " << *eigen_system.second[i] <<";\n";
	}
	return 0;
}
