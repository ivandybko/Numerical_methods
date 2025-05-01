#include <iostream>
#include <fstream>
#include "poisson_equation_2d.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/TridiagonalMatrixAlgorithm.h"

void exportData(const std::vector<std::vector<std::vector<double>>>& data, int pass, const std::string& filename, double tau, double h1, double h2, std::pair<double, double> upper_right_corner, double t0) {
	if (data.empty()) return;

	std::ofstream outFile(filename);
	if (!outFile) {
		std::cerr << "Error opening file!" << '\n';
		return;
	}

	size_t num_points_in_space_1 = std::round(upper_right_corner.first / h1)+1;
	size_t num_points_in_space_2 = std::round(upper_right_corner.second / h2)+1;
	size_t num_points_over_time = std::round(t0 / tau)+1;
	int k = 0;
	for (size_t time = 0; time < num_points_over_time; ++time) {
		for (size_t point = 0; point < num_points_in_space_1; ++point)
		{
			for (size_t point2 = 0; point2 < num_points_in_space_2; ++point2)
			{
				outFile << time*tau << " " << point*h1 << " " << point2*h2 << " " << data[time][point][point2] << '\n';
			}
		}
		time+=pass;
	}

	outFile.close();
}

std::function<double(double,double)> const f_ex1 = [](double x1, double x2){return 0.0;};
BoundaryCondition<double> bc_ex1 = {[](double x){return 1.0;}, false};
auto f_ex2 = f_ex1;
BoundaryCondition<double> bottom_bc_ex2 = {[](double x){return -1.0;}, true};
BoundaryCondition<double> top_bc_ex2 = {[](double x){return 1.0;}, true};
BoundaryCondition<double> left_bc_ex2 = {[](double x){return 1.0+x;}, false};
BoundaryCondition<double> right_bc_ex2 = {[](double x){return 1.0+x;}, false};

std::function<double(double,double)> const f_ex3 = [](double x1, double x2){return 4.0;};
BoundaryCondition<double> left_bc_ex3 = {[](double x){return 0.0;}, true};
BoundaryCondition<double> right_bc_ex3 = {[](double x){return 2.0;}, true};
BoundaryCondition<double> bottom_bc_ex3 = {[](double x){return x*x;}, false};
BoundaryCondition<double> top_bc_ex3 = {[](double x){return 1.0+x*x;}, false};

auto f_ex4 = f_ex1;
BoundaryCondition<double> left_bc_ex4 = {[](double x){return -1.0;}, true};
BoundaryCondition<double> right_bc_ex4 = {[](double x){return 1.0;}, true};
BoundaryCondition<double> bottom_bc_ex4 = {[](double x){return 1.0+x;}, false};
BoundaryCondition<double> top_bc_ex4 = {[](double x){return 1.0+x;}, false};

std::function<double(double,double)> const f_t1 = [](double x1, double x2){return 0.0;};
BoundaryCondition<double> bottom_bc_t1 = {[](double x){return 0.0;}, false};
BoundaryCondition<double> top_bc_t1 = {[](double x){return 1.0;}, false};
BoundaryCondition<double> left_bc_t1 = {[](double x){return x;}, false};
BoundaryCondition<double> right_bc_t1 = {[](double x){return x;}, false};


std::function<double(double,double)> const f_ex5 = [](double x1, double x2){return 4.0;};
BoundaryCondition<double> bottom_bc_ex5 = {[](double x){return 0.0;}, true};
BoundaryCondition<double> top_bc_ex5 = {[](double x){return 2.0;}, true};
BoundaryCondition<double> left_bc_ex5 = {[](double x){return x*x;}, false};
BoundaryCondition<double> right_bc_ex5 = {[](double x){return 1.0+x*x;}, false};

int main()
{
	double tau=0.5;
	double h1=0.1;
	double h2=0.1;
	double t0{1000.0};
	std::pair<double,double> M={1.0,1.0};
	// auto sol1 = poisson_equation_2d<double>(M, t0, f_ex1, bc_ex1, bc_ex1, bc_ex1, bc_ex1, tau, h1, h2);
	// std::string path1 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab4 (Numerical solution of boundary value problems for the two-dimensional Poisson equation)/data/solution1.txt";
	// exportData(sol1, 1000, path1, tau, h1, h2, M ,t0);
	//
	//
	auto sol2 = poisson_equation_2d<double>(M, t0, f_ex2, left_bc_ex2, right_bc_ex2, bottom_bc_ex2, top_bc_ex2, tau, h1, h2);
	std::string path2 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab4 (Numerical solution of boundary value problems for the two-dimensional Poisson equation)/data/solution2.txt";
	exportData(sol2, 10, path2, tau, h1, h2, M ,t0);

	// auto sol3 = poisson_equation_2d<double>(M, t0, f_ex3, left_bc_ex3, right_bc_ex3, bottom_bc_ex3, top_bc_ex3, tau, h1, h2);
	// std::string path3 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab4 (Numerical solution of boundary value problems for the two-dimensional Poisson equation)/data/solution3.txt";
	// exportData(sol3, 100, path3, tau, h1, h2, M ,t0);
	// //
	// //
	// auto solt1 = poisson_equation_2d<double>(M, t0, f_t1, left_bc_t1, right_bc_t1, bottom_bc_t1, top_bc_t1, tau, h1, h2);
	// std::string patht1 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab4 (Numerical solution of boundary value problems for the two-dimensional Poisson equation)/data/t1.txt";
	// exportData(solt1, 100, patht1, tau, h1, h2, M ,t0);

	// auto sol4 = poisson_equation_2d<double>(M, t0, f_ex2, left_bc_ex4, right_bc_ex4, bottom_bc_ex4, top_bc_ex4, tau, h1, h2);
	// std::string path4 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab4 (Numerical solution of boundary value problems for the two-dimensional Poisson equation)/data/solution4.txt";
	// exportData(sol4, 1000, path4, tau, h1, h2, M ,t0);

	// auto sol5 = poisson_equation_2d<double>(M, t0, f_ex5, left_bc_ex5, right_bc_ex5, bottom_bc_ex5, top_bc_ex5, tau, h1, h2);
	// std::string path5 = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab4 (Numerical solution of boundary value problems for the two-dimensional Poisson equation)/data/solution5.txt";
	// exportData(sol5, 1000, path5, tau, h1, h2, M ,t0);
	// return 0;
}

// TIP See CLion help at <a
// href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>.
//  Also, you can try interactive lessons for CLion by selecting
//  'Help | Learn IDE Features' from the main menu.