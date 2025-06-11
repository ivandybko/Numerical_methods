#ifndef QUADRATURE_FREDHOLM_SOLVER_H
#define QUADRATURE_FREDHOLM_SOLVER_H
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/Gauss.h"

#include <vector>
#include <functional>
template <typename T>
std::vector<T> quadrature_fredholm_solve(T lambda, T a, T b, std::function<T(T, T)> K, std::function<T(T)> f, int grid_size){
	  std::vector<std::vector<T>> matrix(grid_size, std::vector<T>(grid_size));
      std::vector<T> rhs(grid_size);
      T h = (b - a) / (grid_size - 1);
     for(int i = 0; i < grid_size; i++){
	 	rhs[i]=f(a+i*h);
	 	matrix[i][0] = 1;
     	for(int k = 0; k < grid_size; k++){
	   		T w_j;
	   		if (k == 0 or k == grid_size - 1)
       		{
       			w_j = h / 2;
	   		}
	        else {
       			w_j = h;
		    }
		    matrix[i][k] = - lambda * w_j * K(i*h, k*h);
        }
        matrix[i][i] += 1;
     }
     auto result = Gauss<T>({matrix, rhs});
     return result;
}
#endif //QUADRATURE_FREDHOLM_SOLVER_H
