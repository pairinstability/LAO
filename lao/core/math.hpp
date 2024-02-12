/// math.hpp implements mathematical expressions on matrix and scalar types.
/// This includes:
/// - `+` addition of two matrices.
/// - `-` subtraction of two matrices, or negation.
/// - `*` matrix multiplication.
/// - `%` element-wise multiplication.
/// - `==` element-wise equality evaluation.
/// - `!=` element-wise non-equality evaluation.
/// - `>=` element-wise greater than or equal to evaluation.
/// - `<=` element-wise less than or equal to evaluation.
/// - `>` element-wise greater than evaluation.
/// - `<` element-wise less than evaluation.

#ifndef LAO_CORE_MATH_H_
#define LAO_CORE_MATH_H_

#include <lao/core/expression.hpp>
#include <lao/core/forward.hpp>
#include <iostream>

namespace lao {

/// @brief Matrix addition.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixAddition : public MatrixExpression<MatrixAddition<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixAddition(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) + m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator+ overload for matrix addition.
template <typename S, size_t R, size_t C, typename E1, typename E2>
auto operator+(const MatrixExpression<E1, S, R, C>& lhs, const MatrixExpression<E2, S, R, C>& rhs)
{
    return MatrixAddition<S, R, C, MatrixExpression<E1, S, R, C>, MatrixExpression<E2, S, R, C>>(lhs, rhs);
}

/// @brief Matrix subtraction.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixSubtraction : public MatrixExpression<MatrixSubtraction<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixSubtraction(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) - m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator- overload for matrix subtraction.
template <typename S, size_t R, size_t C, typename E1, typename E2>
auto operator-(const MatrixExpression<E1, S, R, C>& lhs, const MatrixExpression<E2, S, R, C>& rhs)
{
    return MatrixSubtraction<S, R, C, MatrixExpression<E1, S, R, C>, MatrixExpression<E2, S, R, C>>(lhs, rhs);
}

/*
/// @brief concept for verifying two Scalar types are the same.
template <typename T1, typename T2>
concept ValidScalar = std::is_same_v<T1, T2>;

/// @brief concept for verifying two Scalar types are the same, and the matrices have
/// the same shape.
/// @details this constraint is used for matrix addition and subtraction operations.
template <typename MatrixType1, typename MatrixType2>
concept CompatibleMatrixAdd = requires(MatrixType1 m1, MatrixType2 m2)
{
    requires ValidScalar<typename MatrixType1::value_type, typename MatrixType2::value_type>;
    requires MatrixType1::rows()
    == MatrixType2::rows() && MatrixType1::cols() == MatrixType2::cols();
};
*/

}; // namespace lao

#endif // LAO_CORE_MATH_H_