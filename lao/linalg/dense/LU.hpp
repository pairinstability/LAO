#ifndef LAO_LINALG_DENSE_LU_H_
#define LAO_LINALG_DENSE_LU_H_

#include <cstddef>
#include <lao/linalg/dense/matrix.hpp>
#include <lao/linalg/math/constraints.hpp>

namespace lao {
namespace linalg {

    /// @brief LU decomposition with Doolittle's algorithm.
    /// @details This decomposes a coefficient matrix A into the upper and lower matrices
    /// given the equation Ax = b.
    /// @param A A matrix used in the equation Ax=b to be solved.
    /// @param L Lower diagonal matrix, returned.
    /// @param U Upper diagonal matrix, returned.
    template <typename S, size_t R, size_t C>
    requires EnforceSquareMatrix<S, R, C>
    void LU_doolittle(const Matrix<S, R, C>& A, Matrix<S, R, C>& L, Matrix<S, R, C>& U)
    {
        // initialize L as identity matrix and U as zero matrix
        L.eye();
        U.zeros();

        for (size_t j = 1; j < R + 1; ++j) {
            for (size_t i = j; i < C + 1; ++i) {
                S sum = 0;
                for (size_t k = 1; k < j; ++k)
                    sum += L(j, k) * U(k, i);
                U(j, i) = A(j, i) - sum;
            }

            for (size_t i = j + 1; i < R + 1; ++i) {
                S sum = 0;
                for (size_t k = 1; k < j; ++k)
                    sum += L(i, k) * U(k, j);
                L(i, j) = (A(i, j) - sum) / U(j, j);
            }
        }
    };

}; // namespace linalg
}; // namespace lao

#endif // LAO_LINALG_DENSE_LU_H