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

#include <iostream>
#include <lao/core/expression.hpp>
#include <lao/core/forward.hpp>

namespace lao {

/// @brief Concept for enforcing matrix addition/subtraction constraints.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
concept MatrixCompatible = requires
{
    // Check if the scalar types match
    requires std::same_as<S1, S2>;
    // Check if the dimensions match
    requires R1 == R2&& C1 == C2;
};

/// @brief Concept for enforcing matrix multiplication constraints.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
concept MatrixCompatibleMult = requires
{
    // Check if the scalar types match
    requires std::same_as<S1, S2>;
    // Check if the dimensions match
    requires C1 == R2;
};

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
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires MatrixCompatible<S1, S2, R1, R2, C1, C2, E1, E2>
auto operator+(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixAddition<S1, R1, C1, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
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
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
auto operator-(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixSubtraction<S1, R1, C1, MatrixExpression<E1, S2, R2, C2>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix multiplication.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixMultiplication : public MatrixExpression<MatrixMultiplication<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixMultiplication(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        S dot_product = S();
        for (size_t i = 0; i < m_lhs.cols(); ++i) {
            dot_product += m_lhs(row, i) * m_rhs(i, col);
        }
        return dot_product;
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator* overload for matrix multiplication.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires MatrixCompatibleMult<S1, S2, R1, R2, C1, C2, E1, E2>
auto operator*(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixMultiplication<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

}; // namespace lao

#endif // LAO_CORE_MATH_H_