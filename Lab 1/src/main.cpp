#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include "Gauss.h"
#include "QR decomposition.h"
template <typename T>
std::unique_ptr<std::pair<std::vector<std::vector<T>>, std::vector<T>>> readData(const std::string& path) {
	auto data = std::make_unique<std::pair<std::vector<std::vector<T>>, std::vector<T>>>();
	std::vector<std::vector<T>>& matrix = data->first;
	std::vector<T>& vector = data->second;
	std::ifstream inputMatrix(path + "/matrix.dat");
	if (!inputMatrix.is_open()) {
		std::cerr << "Файл matrix.dat не найден" << std::endl;
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
		std::cerr << "Файл vector.dat не найден" << std::endl;
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
std::vector<T> readSolution(const std::string& path) {
	std::ifstream inputSolution(path + "/solution.dat");
	if (!inputSolution.is_open()) {
		std::cerr << "Файл solution.dat не найден. Проверка решения невозможна" << std::endl;
		return {0};
	}
	int n;
	inputSolution >> n;
	std::vector<T> solution(n);
	for (auto& element : solution) {
		inputSolution >> element;
	}
	inputSolution.close();
	return solution;
}

int main(){
	std::cout << "Введите путь к директории: ";
	std::string path;
	std::cin >> path;
	auto data = *readData<float>(path);
	auto x_gauss = Gauss(data);
	auto solution= readSolution<float>(path);
	std::cout << "Решение методом Гаусса:\n { " ;
	for (auto element : x_gauss){
		std::cout << element << " ";
	}
	std::cout << "}\n";
	if (solution[0] !=0){
		if (x_gauss == solution){
			std::cout << "Решение методом Гаусса найдено верно";
		}
		else{
			std::cerr << "Решение методом Гаусса найдено неверно";
			std::cout << "Истинное решение:\n { " ;
			for (auto element : solution){
				std::cout << element <<" ";
			}
			std::cout << "}\n";
		}
	}
	return 0;
}