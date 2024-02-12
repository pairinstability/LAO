/// matrix.hpp implements the Matrix class as a class for dense matrix operations, and
/// the partial specializations of this to implement row and column vectors.

#ifndef LAO_DENSE_MATRIX_H_
#define LAO_DENSE_MATRIX_H_

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <lao/core/expression.hpp>
#include <lao/core/forward.hpp>
#include <lao/core/math.hpp>
#include <random>
#include <stdexcept>
#include <vector>

namespace lao {

enum class filltype {
    zeros,
    ones,
    eye,
    rand,
    none
};

/// @brief The class representing a dense matrix.
/// @details A matrix is represented by a number of rows and columns,
/// noted as `row x column`.
/// It is templated with a Scalar parameter, and a Row and Column.
/// The default storage mechanism is using std::vector with row-major ordering,
/// though this can be swapped out assuming it has a linear access pattern.
template <typename S, size_t R, size_t C, typename B>
class Matrix : public MatrixExpression<Matrix<S, R, C>, S, R, C> {
public:
    using value_type = S;
    using storage_type = B;

    /// @brief Default constructor which zero initializes the matrix.
    Matrix()
        : m_elements(R * C, value_type(0))
    {
    }

    /// @brief Constructor taking a vector of elements to construct a matrix.
    Matrix(const std::vector<S>& elements)
    {
        if (elements.size() != R * C)
            throw std::invalid_argument("Initializer list does not match matrix size.");
        m_elements = elements;
    }

    /// @brief Constructor to initialize with an initializer list.
    Matrix(std::initializer_list<std::initializer_list<value_type>> list)
        : m_elements(R * C, value_type(0))
    {
        if (list.size() != R || list.begin()->size() != C)
            throw std::invalid_argument("Initializer list does not match matrix size.");
        size_t row = 0;
        for (const auto& il : list) {
            std::copy(il.begin(), il.end(), m_elements.begin() + row * C);
            ++row;
        }
    }

    /// @brief Constructor to initialize the matrix with a fill type.
    Matrix(filltype fill)
        : m_elements(R * C, value_type(0))
    {
        if (fill == filltype::zeros) {
            // Set the matrix to all zeros
            zeros();
        } else if (fill == filltype::ones) {
            // Set the matrix to all ones
            ones();
        } else if (fill == filltype::eye) {
            // Set the matrix to the identity, i.e. the diagonal has ones.
            // This is only valid for square matrics.
            eye();
        } else if (fill == filltype::rand) {
            // Fill in the matrix with random values in range [0,1].
            rand();
        } else if (fill == filltype::none) {
            // nada
        } else {
            throw std::invalid_argument("Invalid filltype.");
        }
    }

    /// @brief Copy constructor.
    Matrix(const Matrix& other)
        : m_elements(other.m_elements)
    {
    }

    Matrix& operator=(const Matrix& other)
    {
        if (&other != this)
            m_elements = other.m_elements;
        return *this;
    }

    /// @brief Operator for converting MatrixExpression <-> Matrix
    template <typename E>
    Matrix(const MatrixExpression<E, S, R, C>& expr)
        : m_elements(R * C, value_type(0))
    {
        for (size_t i = 0; i < R; ++i)
            for (size_t j = 0; j < C; ++j)
                m_elements[i * C + j] = static_cast<value_type>(expr(i, j));
    }

    /// @brief operator overload for () to access elements.
    value_type& operator()(size_t row, size_t col)
    {
        if (row > R || col > C)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row * C + col];
    }

    /// @brief operator overload for () to access elements.
    const value_type& operator()(size_t row, size_t col) const
    {
        if (row > R || col > C)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row * C + col];
    }

    /// @brief Returns the number of rows.
    size_t rows() const noexcept
    {
        return R;
    }

    /// @brief Returns the number of columns.
    size_t cols() const noexcept
    {
        return C;
    }

    /// @brief Checks if the matrix is empty.
    bool is_empty() const
    {
        return m_elements.empty();
    }

    /// @brief Sets all elements to zero.
    void zeros()
    {
        std::fill(m_elements.begin(), m_elements.end(), value_type(0));
    }

    /// @brief Sets all elements to ones.
    void ones()
    {
        std::fill(m_elements.begin(), m_elements.end(), value_type(1));
    }

    /// @brief Sets all elements to the identity matrix.
    void eye()
    {
        if constexpr (R == C) {
            std::fill(m_elements.begin(), m_elements.end(), value_type(0));
            for (size_t i = 0; i < R; ++i) {
                m_elements[i * C + i] = value_type(1);
            }
        } else {
            throw std::logic_error("Identity matrix is only defined for square matrices.");
        }
    }

    /// @brief Sets all elements to random values in the [0, 1] interval.
    void rand()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<value_type> dis(0.0, 1.0);
        for (size_t i = 0; i < m_elements.size(); ++i)
            m_elements[i] = dis(gen);
    }

    /// @brief Sets all elements to a specified value.
    void fill(value_type val)
    {
        std::fill(m_elements.begin(), m_elements.end(), val);
    }

    /// @brief Sets all elements using a lambda function.
    void fillf(std::function<value_type()> lambda)
    {
        for (size_t i = 0; i < m_elements.size(); ++i)
            m_elements[i] = lambda();
    }

    /// @brief Resets the matrix to empty.
    void reset()
    {
        m_elements.clear();
    }

    /// @brief Printing method implementation, for printing the matrix.
    friend std::ostream& operator<<(std::ostream& os, const Matrix<S, R, C>& matrix)
    {
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < C; ++j) {
                os << matrix(i, j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

private:
    storage_type m_elements;
};

}; // namespace lao

#endif // LAO_DENSE_MATRIX_H_