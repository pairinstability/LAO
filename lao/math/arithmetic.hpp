#ifndef LAO_MATH_ARITHMETIC_H_
#define LAO_MATH_ARITHMETIC_H_

#include <lao/core/expression.hpp>
#include <lao/math/constraints.hpp>

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
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
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
requires EnforceMatMulReqs<S1, S2, R1, R2, C1, C2>
auto operator*(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixMultiplication<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Scalar-Matrix multiplication.
template <typename S, size_t R, size_t C, typename E>
class MatrixScalarMultiplication : public MatrixExpression<MatrixScalarMultiplication<S, R, C, E>, S, R, C> {
public:
    using left_expr = S;
    using right_expr = E;

    MatrixScalarMultiplication(const S& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs * m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator* overload for scalar-matrix multiplication.
template <typename S1, typename S2, size_t R, size_t C, typename E>
requires EnforceSameType<S1, S2>
auto operator*(const MatrixExpression<E, S1, R, C>& lhs, const S2& rhs)
{
    return MatrixScalarMultiplication<S2, R, C, MatrixExpression<E, S1, R, C>>(lhs, rhs);
}

template <typename S1, typename S2, size_t R, size_t C, typename E>
requires EnforceSameType<S1, S2>
auto operator*(const S1& lhs, const MatrixExpression<E, S2, R, C>& rhs)
{
    return MatrixScalarMultiplication<S1, R, C, MatrixExpression<E, S2, R, C>>(lhs, rhs);
}

/// @brief Matrix element-wise multiplication.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseMultiplication : public MatrixExpression<MatrixElementWiseMultiplication<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseMultiplication(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) * m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator% overload for matrix element-wise multiplication.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator%(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseMultiplication<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix element-wise equality check.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseEquality : public MatrixExpression<MatrixElementWiseEquality<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseEquality(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) == m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator== overload for matrix element-wise equality check.
/// @details If the elements are the same, element in return matrix is 1. Else, its 0.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator==(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseEquality<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix element-wise non-equality check.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseNonEquality : public MatrixExpression<MatrixElementWiseNonEquality<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseNonEquality(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) != m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator!= overload for matrix element-wise non equality check.
/// @details If the elements are the same, element in return matrix is 0. Else, its 1.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator!=(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseNonEquality<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix element-wise greater than equal check.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseGEQ : public MatrixExpression<MatrixElementWiseGEQ<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseGEQ(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) >= m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator>= overload for matrix element-wise greater than equal check.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator>=(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseGEQ<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix element-wise greater than check.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseGT : public MatrixExpression<MatrixElementWiseGT<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseGT(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) > m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator> overload for matrix element-wise greater than check.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator>(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseGT<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix element-wise greater less than equal check.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseLEQ : public MatrixExpression<MatrixElementWiseLEQ<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseLEQ(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) <= m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator<= overload for matrix element-wise less than equal check.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator<=(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseLEQ<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

/// @brief Matrix element-wise greater less than check.
template <typename S, size_t R, size_t C, typename E1, typename E2>
class MatrixElementWiseLT : public MatrixExpression<MatrixElementWiseLT<S, R, C, E1, E2>, S, R, C> {
public:
    using left_expr = E1;
    using right_expr = E2;

    MatrixElementWiseLT(const left_expr& lhs, const right_expr& rhs)
        : m_lhs(lhs)
        , m_rhs(rhs)
    {
    }

    S operator()(size_t row, size_t col) const
    {
        return m_lhs(row, col) <= m_rhs(row, col);
    }

private:
    const left_expr& m_lhs;
    const right_expr& m_rhs;
};

/// @brief operator<  overload for matrix element-wise less than check.
template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2, typename E1, typename E2>
requires EnforceSameShape<S1, S2, R1, R2, C1, C2>
auto operator<(const MatrixExpression<E1, S1, R1, C1>& lhs, const MatrixExpression<E2, S2, R2, C2>& rhs)
{
    return MatrixElementWiseLT<S1, R1, C2, MatrixExpression<E1, S1, R1, C1>, MatrixExpression<E2, S2, R2, C2>>(lhs, rhs);
}

};

#endif // LAO_MATH_ARITHMETIC_H_