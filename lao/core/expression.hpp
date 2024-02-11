#ifndef LAO_CORE_BASE_H_
#define LAO_CORE_BASE_H_

#include <cstddef>
#include <iostream>

namespace lao {

/// @brief MatrixExpression represents a matrix expression.
/// @details This base class represents an expression when building expression
/// trees for matrix mathematics operations. Conceptually, this is an expression
/// template used for deferred evaluation of expressions which allows for efficient
/// expression composition without immediate computation of intermediate results.
template <typename Scalar, typename Derived>
class MatrixExpression {
public:
    /// @brief Returns the number of rows from the derived class.
    size_t rows() const
    {
        return static_cast<const Derived*>(this)->rows();
    }

    /// @brief Returns the number of cols from the derived class.
    size_t cols() const
    {
        return static_cast<const Derived*>(this)->cols();
    }

    /// @brief Returns the shape of the matrix.
    size_t shape() const noexcept
    {
        return rows() * cols();
    }

    /// @brief Prints the contents of the matrix to stdout.
    /// @details Caveat emptor.
    void print(std::ostream& os) const
    {
        static_cast<const Derived*>(this)->print(os);
    }

    /// @brief operator overload for () to access elements.
    Scalar operator()(size_t row, size_t col) const
    {
        return static_cast<const Derived*>(this)->operator()(row, col);
    }
};

template <typename Scalar, typename Derived>
std::ostream& operator<<(std::ostream& os, const MatrixExpression<Scalar, Derived>& derived)
{
    derived.print(os);
    return os;
}

}; // namespace lao

#endif // LAO_CORE_BASE_H_