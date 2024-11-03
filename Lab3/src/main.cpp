#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "FrancisKublanovskayaQR.h"
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
	auto eigen_vectors_QR = QRAlgorithm<double>(*data, 1e-4, true);
	std::cout << eigen_vectors_QR;
	return 0;
}
