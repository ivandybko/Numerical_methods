#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <algorithm>
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
	std::cout << "Решение во float:\n";
	auto data_f = *readData<float>(path);
	auto x_gauss_f = Gauss<float>(data_f);
	auto solution_f= readSolution<float>(path);
	if (!x_gauss_f.empty()){
		std::cout << "Решение методом Гаусса:\n { " ;
		for (auto element : x_gauss_f){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		if (x_gauss_f == solution_f){
			std::cout << "Решение методом Гаусса найдено верно";
		}
		else{
			std::cerr << "Решение методом Гаусса найдено неверно";
			std::cout << "Истинное решение:\n { " ;
			for (auto element : solution_f){
				std::cout << element <<" ";
			}
			std::cout << "}\n";
		}
	}
	else{
		std::cerr << "Решения не существует или оно не единственно";
	}
	auto x_QR = QRdecompostition(data_f);
	std::cout << "Решение в double:\n";
	auto data_d = *readData<double>(path);
	auto x_gauss_d = Gauss<double>(data_d);
	auto solution_d= readSolution<double>(path);
	if (!x_gauss_d.empty()){
		std::cout << "Решение методом Гаусса:\n { " ;
		for (auto element : x_gauss_d){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		if (x_gauss_d == solution_d){
			std::cout << "Решение методом Гаусса найдено верно";
		}
		else{
			std::cerr << "Решение методом Гаусса найдено неверно";
			std::cout << "Истинное решение:\n { " ;
			for (auto element : solution_d){
				std::cout << element <<" ";
			}
			std::cout << "}\n";
		}
	}
	else{
		std::cerr << "Решения не существует или оно не единственно";
	}
	auto x_QR_d = QRdecompostition(data_d);
	return 0;
}