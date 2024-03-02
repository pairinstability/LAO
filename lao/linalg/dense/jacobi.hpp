#ifndef LAO_LINALG_DENSE_JACOBI_H_
#define LAO_LINALG_DENSE_JACOBI_H_

#include <lao/linalg/dense/matrix.hpp>
#include <lao/linalg/math/constraints.hpp>

namespace lao {
namespace linalg {

    /// @brief The Jacobi method for solving a system of linear equations.
    /// @details Given Ax = b where A is a known square matrix, b is a known column
    /// matrix, and x must be found. The Jacobi method can be used to approximate x.
    /// This uses the element-based formula.
    /// @param x The n x 1 matrix to solve for, x.
    /// @param A Known n x n matrix, A.
    /// @param b Known n x 1 matrix, b.
    /// @param max_iterations Maximum number of iterations.
    /// @param tol The numerical tolerance.
    template <typename S, size_t R, size_t C>
    requires EnforceSquareMatrix<S, R, C>
    void solve_jacobi_element(Matrix<S, R, 1>& x, const Matrix<S, R, C>& A, Matrix<S, R, 1>& b, size_t max_iterations, S tol)
    {
        // initialize solution vector to initial guess. Matrix default
        // constructor initializes matrix to all zeros.
        x.ones();
        Matrix<S, R, 1> x_scratch;

        for (size_t iter = 0; iter < max_iterations; ++iter) {
            S err = S(0);

            for (size_t i = 1; i < R + 1; ++i) {
                S sum = S(0);
                for (size_t j = 1; j < R + 1; ++j) {
                    if (j != i)
                        sum += A(i, j) * x(j, 1);
                }

                x_scratch(i, 1) = (b(i, 1) - sum) / A(i, i);
            }

            for (size_t i = 1; i < R + 1; ++i) {
                err += std::abs(x_scratch(i, 1) - x(i, 1));
            }

            if (err < tol) {
                std::cout << "conv" << std::endl;
                x = x_scratch;
                return;
            }

            x = x_scratch;
        }

        std::cout << "no conv" << std::endl;
        return;
    };

}; // namespace linalg
}; // namespace lao

#endif // LAO_LINALG_DENSE_JACOBI_H_