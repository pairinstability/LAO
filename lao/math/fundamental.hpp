#ifndef LAO_MATH_FUNDAMENTAL_H_
#define LAO_MATH_FUNDAMENTAL_H_

#include <lao/core/expression.hpp>
#include <lao/math/constraints.hpp>

namespace lao {

/// @brief Accumulates all of the elements in a matrix.
template <typename S, size_t R, size_t C, typename E>
S accumulate(const MatrixExpression<E, S, R, C>& matrix)
{
    S sum = 0;
    for (size_t i = 1; i < matrix.rows()+1; ++i) {
        for (size_t j = 1; j < matrix.cols()+1; ++j) {
            sum += matrix(i, j);
        }
    }
    return sum;
}

/// @brief Calculates the determinant of a matrix.
template <typename S, size_t R, size_t C, typename E>
S det(const MatrixExpression<E, S, R, C>& matrix)
{
    // TODO
}

/// @brief Calculates the rank of a matrix.
template <typename S, size_t R, size_t C, typename E>
size_t rank(const MatrixExpression<E, S, R, C>& matrix)
{
    // TODO
}

template <typename S, size_t R, size_t C, typename E>
class MatrixTranspose : public MatrixExpression<MatrixTranspose<S, R, C, E>, S, C, R> {
public:
    using value_type = S;

    MatrixTranspose(const MatrixExpression<E, S, R, C>& matrix)
        : m_matrix(matrix)
    {
    }

    S operator()(size_t row, size_t col) const override
    {
        return m_matrix(col, row);
    }

private:
    const MatrixExpression<E, S, R, C>& m_matrix;
};

/// @brief Returns the transposition of a matrix.
template <typename S, size_t R, size_t C, typename E>
auto transpose(const MatrixExpression<E, S, R, C>& matrix)
{
    return MatrixTranspose<S, R, C, E>(matrix);
}

/// @brief Calculates the trace of a matrix.
/// @details Only works for square matrices.
template <typename S, size_t R, size_t C, typename E>
requires EnforceSquareMatrix<S, R, C>
    S trace(const MatrixExpression<E, S, R, C>& matrix)
{
    S sum = 0;
    for (size_t i = 1; i < std::min(matrix.rows()+1, matrix.cols()+1); ++i) {
        sum += matrix(i, i);
    }
    return sum;
}

template <typename S, size_t R, size_t C, typename E>
class MatrixInverse : public MatrixExpression<MatrixInverse<S, R, C, E>, S, C, R> {
public:
    using value_type = S;

    MatrixInverse(const MatrixExpression<E, S, R, C>& matrix)
        : m_matrix(matrix)
    {
    }

    S operator()(size_t row, size_t col) const override
    {

    }

private:
    const MatrixExpression<E, S, R, C>& m_matrix;
};

/// @brief Calculates the inverse of a matrix.
/// @details Only works for square matrices.
template <typename S, size_t R, size_t C, typename E>
requires EnforceSquareMatrix<S, R, C>
    auto inv(const MatrixExpression<E, S, R, C>& matrix)
{
    // TODO
}

};

#endif // LAO_MATH_FUNDAMENTAL_H_