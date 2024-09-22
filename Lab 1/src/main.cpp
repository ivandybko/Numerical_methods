#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
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
		std::cerr << "Файл solution.dat не найден. Проверка решения невозможна\n" << std::endl;
		return std::vector<T>();
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
template <typename T>
std::pair<T,T> residual(const std::vector<std::vector<T>>& A, const std::vector<T> & b_true, const std::vector<T> &x){
	size_t n = A.size();
	std::vector<T> b(n,0);
	T norm1{0},normi{0};
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			b[i]+=A[i][j]*x[j];
		}
		norm1+=abs(b_true[i]-b[i]);
		if (abs(b_true[i]-b[i]) > normi){
			normi=abs(b_true[i]-b[i]);
		}
	}
	return std::pair<T,T>(norm1,normi);
}

int main(){
	std::cout << "Введите путь к директории: ";
	std::string path;
	std::cin >> path;
	std::cout << "\nРешение во float:\n";
	auto data_f = *readData<float>(path);
	auto x_gauss_f = Gauss<float>(data_f);
	auto solution_f= readSolution<float>(path);
	if (!solution_f.empty()){
		std::cout << "Точное решение:\n { " ;
		for (auto element : solution_f){
			std::cout << element <<" ";
		}
		std::cout << "}\n";
	}
	if (!x_gauss_f.empty()){
		std::cout << "Решение методом Гаусса:\n { " ;
		for (auto element : x_gauss_f){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_f = residual<float>(data_f.first, data_f.second, x_gauss_f);
		std::cout << "Октаэдрическая норма: " << norm_f.first << ';' << " Кубическая норма: " << norm_f.second << '\n';
	}
	else{
		std::cout << "Метод Гаусса: Решения не существует или оно не единственно\n";
	}

	auto x_QR_f = QRdecompostition<float>(data_f);
	if (!x_QR_f.empty()){
		std::cout << "Решение используя QR разложение:\n { " ;
		for (auto element : x_QR_f){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_f = residual<float>(data_f.first, data_f.second, x_QR_f);
		std::cout << "Октаэдрическая норма: " << norm_f.first << ';' << " Кубическая норма: " << norm_f.second << '\n';
	}
	else{
		std::cout << "Метод QR разложения: Решения не существует или оно не единственно\n";
	}

	std::cout << "\nОценка числа обусловленности:\n" << "cond1=" << (x_gauss_f.empty() ? std::numeric_limits<double>::infinity() :octahedralNorm<float>(data_f.first)* octahedralNorm<float>(
		*inverse<float>(data_f.first))) << '\t' << "cond∞=" << (x_gauss_f.empty() ? std::numeric_limits<double>::infinity() : cubicNorm<float>(data_f.first)* cubicNorm<float>(
		*inverse<float>(data_f.first))) << '\n';


	std::cout << "\nРешение в double:\n";
	auto data_d = *readData<double>(path);
	auto x_gauss_d = Gauss<double>(data_d);
	auto solution_d= readSolution<double>(path);
	if (!x_gauss_d.empty()){
		std::cout << "Решение методом Гаусса:\n { " ;
		for (auto element : x_gauss_d){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_d = residual<double>(data_d.first, data_d.second, x_gauss_d);
		std::cout << "Октаэдрическая норма: " << norm_d.first << ';' << " Кубическая норма: " << norm_d.second << '\n';
	}
	else{
		std::cout << "Метод Гаусса: Решения не существует или оно не единственно\n";
	}

	auto x_QR_d = QRdecompostition<double>(data_d);
	if (!x_QR_d.empty()){
		std::cout << "Решение используя QR разложение:\n { " ;
		for (auto element : x_QR_d){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_d = residual<double>(data_d.first, data_d.second, x_QR_d);
		std::cout << "Октаэдрическая норма: " << norm_d.first << ';' << " Кубическая норма: " << norm_d.second << '\n';
	}
	else{
		std::cout << "Метод QR разложения: Решения не существует или оно не единственно\n";
	}

	std::cout << "\nОценка числа обусловленности:\n" << "cond1=" << (x_gauss_d.empty() ? std::numeric_limits<double>::infinity() : octahedralNorm<double>(data_d.first)* octahedralNorm<double>(
		*inverse<double>(data_d.first))) << '\t' << "cond∞=" << (x_gauss_d.empty() ? std::numeric_limits<double>::infinity() : cubicNorm<double>(data_d.first)* cubicNorm<double>(
		*inverse<double>(data_d.first))) << '\n';

	size_t n{data_f.first.size()};
	std::vector<float> b_f_per(data_f.second);
	for (int i = 0; i < n/2; ++i)
	{
		b_f_per[rand() % n] += 0.01;
		b_f_per[rand() % n] -= 0.01 ;
	}
	std::cout << "\nВектор с внесенным возмущением (в float):\n{ ";
	for (auto element : b_f_per)
	{
		std::cout << element << " ";
	}
	std::cout << " }\n";

	auto x_gauss_f_per = Gauss<float>(std::pair<std::vector<std::vector<float>>, std::vector<float>>(data_f.first, b_f_per));
	if (!x_gauss_f_per.empty()){
		std::cout << "Решение методом Гаусса системы с возмущенной правой частью:\n { " ;
		for (auto element : x_gauss_f_per){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_f = residual<float>(data_f.first, b_f_per, x_gauss_f_per);
		std::cout << "Октаэдрическая норма: " << norm_f.first << ';' << " Кубическая норма: " << norm_f.second << '\n';
	}
	else{
		std::cout << "Метод Гаусса: Решения не существует или оно не единственно\n";
	}

	auto x_QR_f_per = QRdecompostition<float>(std::pair<std::vector<std::vector<float>>, std::vector<float>>(data_f.first, b_f_per));
	if (!x_QR_f_per.empty()){
		std::cout << "Решение используя QR разложение с возмущенной правой частью:\n { " ;
		for (auto element : x_QR_f_per){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_f = residual<float>(data_f.first, b_f_per, x_QR_f_per);
		std::cout << "Октаэдрическая норма: " << norm_f.first << ';' << " Кубическая норма: " << norm_f.second << '\n';
	}
	else{
		std::cout << "Метод QR разложения: Решения не существует или оно не единственно\n";
	}

	std::vector<double> b_d_per(data_d.second);
	for (int i = 0; i < n/2; ++i)
	{
		b_d_per[rand() % n] += 0.01;
		b_d_per[rand() % n] -= 0.01 ;
	}
	std::cout << "\nВектор с внесенным возмущением (в double) [могут быть изменены другие элементы]:\n{ ";
	for (auto element : b_d_per)
	{
		std::cout << element << " ";
	}
	std::cout << "}\n";

	auto x_gauss_d_per = Gauss<double>(std::pair<std::vector<std::vector<double>>, std::vector<double>>(data_d.first, b_d_per));
	if (!x_gauss_d_per.empty()){
		std::cout << "Решение методом Гаусса системы с возмущенной правой частью:\n { " ;
		for (auto element : x_gauss_d_per){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_d = residual<double>(data_d.first, b_d_per, x_gauss_d_per);
		std::cout << "Октаэдрическая норма: " << norm_d.first << ';' << " Кубическая норма: " << norm_d.second << '\n';
	}
	else{
		std::cout << "Метод Гаусса: Решения не существует или оно не единственно\n";
	}

	auto x_QR_d_per = QRdecompostition<double>(std::pair<std::vector<std::vector<double>>, std::vector<double>>(data_d.first, b_d_per));
	if (!x_QR_d_per.empty()){
		std::cout << "Решение используя QR разложение с возмущенной правой частью:\n { " ;
		for (auto element : x_QR_d_per){
			std::cout << element << " ";
		}
		std::cout << "}\n";
		auto norm_d = residual<double>(data_d.first, b_d_per, x_QR_d_per);
		std::cout << "Октаэдрическая норма: " << norm_d.first << ';' << " Кубическая норма: " << norm_d.second << '\n';
	}
	else{
		std::cout << "Метод QR разложения: Решения не существует или оно не единственно\n";
	}
	return 0;
}