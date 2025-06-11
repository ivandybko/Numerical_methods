#ifndef SINGULAR_FREDHOLM_SOLVER_H
#define SINGULAR_FREDHOLM_SOLVER_H
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/Gauss.h"
#include "../../../Linear algebra/Lab1 (Direct methods for solving systems of linear equations)/src/matrix_operations.h"
#include <vector>
#include <functional>

template <typename T>
std::vector<T> singular_fredholm_solver(std::function<T(T)> f, int N) {
    std::vector<std::vector<T>> collocation_points(N);
    std::vector<std::vector<T>> quadrature_points(N);
    T delta_l = 2 * M_PI / N;
    for (int i = 1; i <= N; ++i) {
        T theta_c = 2 * M_PI * (i - 0.5) / N;
        T theta_q = 2 * M_PI * (i-1) / N;
        collocation_points[i-1] = {std::cos(theta_c), std::sin(theta_c)};
        quadrature_points[i-1]  = {std::cos(theta_q),  std::sin(theta_q)};
    }
    std::string filename = "/Users/ivandybko/Projects/Numerical_methods/Mathematical physics/Lab5 (Numerical methods for solving integral equations)/data/singular.txt";
    std::ofstream outFile(filename);


    std::vector<std::vector<T>> A(N+1, std::vector<T>(N+1));
    std::vector<T> b(N+1);

    auto Q_kernel = [&](const std::vector<T>& r, const std::vector<T>& rho) -> std::vector<T>{
        auto diff = r - rho;
        auto sq = diff * diff;
        T factor = 1.0 / (2 * M_PI * sq);
        return {-factor * diff[1], factor * diff[0]};
    };

    for (int i = 0; i < N; ++i) {
        auto &r = collocation_points[i];
        auto phi = atan(r[1] / r[0]);
        b[i] = f(phi);
        for (int j = 0; j < N; ++j) {
            const std::vector<T> rho = quadrature_points[j];
            std::vector<T> Q = Q_kernel(r, rho);
            A[i][j] = delta_l * (Q * r);
        }
        A[i][N] = 1;
    }
    for (int j = 0; j < N; ++j) {
        A[N][j] = delta_l;
    }
    b[N] = 0;
    return Gauss<T>({A, b});
}

#endif //SINGULAR_FREDHOLM_SOLVER_H
