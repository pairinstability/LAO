/// matrix.hpp implements the Matrix class as a class for dense matrix operations, and
/// the partial specializations of this to implement row and column vectors. They all
/// derive from a MatrixExpression base class, which represents a fundamental expression
/// so these objects can be used with arithmetic operations.

#ifndef LAO_DENSE_MATRIX_H_
#define LAO_DENSE_MATRIX_H_

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <lao/core/expression.hpp>
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
/// @details This derives from the MatrixExpression class, and is a leaf in the
/// expression tree.
/// A matrix is represented by a number of rows and columns, noted as `row x column`.
/// The storage mechanism is using std::vector with row-major ordering.
template <typename Scalar, size_t Rows, size_t Cols>
class Matrix : public MatrixExpression<Scalar, Matrix<Scalar, Rows, Cols>> {
public:
    using value_type = Scalar;
    using storage_type = std::vector<value_type>;

    /// @brief Default constructor which zero initializes the matrix.
    Matrix()
        : m_elements(Rows * Cols, value_type(0))
    {
    }

    /// @brief Constructor taking a vector of elements to construct a matrix.
    Matrix(const std::vector<Scalar>& elements)
    {
        if (elements.size() != Rows * Cols) throw std::invalid_argument("Initializer list does not match matrix size.");
        m_elements = elements;
    }

    /// @brief Constructor to initialize with an initializer list.
    Matrix(std::initializer_list<std::initializer_list<value_type>> list)
        : m_elements(Rows * Cols, value_type(0))
    {
        if (list.size() != Rows || list.begin()->size() != Cols) throw std::invalid_argument("Initializer list does not match matrix size.");
        size_t row = 0;
        for (const auto& il : list) {
            std::copy(il.begin(), il.end(), m_elements.begin() + row * Cols);
            ++row;
        }
    }

    /// @brief Constructor to initialize the matrix with a fill type.
    Matrix(filltype fill)
        : m_elements(Rows * Cols, value_type(0))
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

    /// @brief Copy assignment operator.
    Matrix& operator=(const Matrix& other)
    {
        if (&other != this)
            m_elements = other.m_elements;
        return *this;
    }

    /// @brief operator overload for () to access elements.
    value_type& operator()(size_t row, size_t col)
    {
        if (row > Rows || col > Cols)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row * Cols + col];
    }

    /// @brief operator overload for () to access elements.
    const value_type& operator()(size_t row, size_t col) const
    {
        if (row > Rows || col > Cols)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row * Cols + col];
    }

    /// @brief Returns the number of rows.
    size_t rows() const noexcept
    {
        return Rows;
    }

    /// @brief Returns the number of columns.
    size_t cols() const noexcept
    {
        return Cols;
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
        if constexpr (Rows == Cols) {
            std::fill(m_elements.begin(), m_elements.end(), value_type(0));
            for (size_t i = 0; i < Rows; ++i) {
                m_elements[i * Cols + i] = value_type(1);
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
    void print(std::ostream& os) const
    {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                os << m_elements[i * Cols + j] << " ";
            }
            os << "\n";
        }
    }

private:
    storage_type m_elements;
};

/// @brief Unique specialization for a 1x1 matrix.
/// @details This would be rare, but the problem without having this specialization is that
/// its ambiguous, given the vectors are specialized on either 1 row or 1 column, and a 1x1
/// matrix is both. Although this is just a 1x1 matrix and it doesn't need to be a vector,
/// its easier to store it like that for math operations although it isn't efficient.
template <typename Scalar>
class Matrix<Scalar, 1, 1> : public MatrixExpression<Scalar, Matrix<Scalar, 1, 1>> {
public:
    using value_type = Scalar;
    using storage_type = std::vector<value_type>;

    /// @brief Default constructor which zero initializes the vector.
    Matrix()
        : m_elements(1, value_type(0))
    {
    }

    /// @brief Constructor taking a vector of elements to construct a vector.
    Matrix(const std::vector<Scalar>& elements)
    {
        static_assert(elements.size() == 1, "Initializer list size does not match matrix size.");
        m_elements = elements;
    }

    /// @brief Constructor to initialize the matrix with a fill type.
    /// @details As this is a vector, it does not support the eye fill.
    Matrix(filltype fill)
        : m_elements(1, value_type(0))
    {
        if (fill == filltype::zeros) {
            // Set the matrix to all zeros
            m_elements[0] = value_type(0);
        } else if (fill == filltype::ones) {
            // Set the matrix to all ones
            m_elements[0] = value_type(1);
        } else if (fill == filltype::rand) {
            // Fill in the matrix with random values in range [0,1].
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<value_type> dis(0.0, 1.0);
            m_elements[0] = dis(gen);
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

    /// @brief Copy assignment operator.
    Matrix& operator=(const Matrix& other)
    {
        if (&other != this)
            m_elements = other.m_elements;
        return *this;
    }

    /// @brief operator overload for () to access elements.
    value_type& operator()(size_t index)
    {
        return m_elements[index];
    }

    /// @brief operator overload for () to access elements.
    const value_type& operator()(size_t index) const
    {
        return m_elements[index];
    }

    /// @brief Operator overload for () to access elements.
    value_type& operator()(size_t row, size_t col)
    {
        if (row != 0 || col != 0)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row];
    }

    /// @brief Operator overload for () to access elements.
    const value_type& operator()(size_t row, size_t col) const
    {
        if (row != 0 || col != 0)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row];
    }

    /// @brief Checks if the vector is empty.
    bool is_empty() const
    {
        return m_elements.empty();
    }

    /// @brief Sets all elements to zero.
    void zeros()
    {
        m_elements[0] = value_type(0);
    }

    /// @brief Sets all elements to ones.
    void ones()
    {
        m_elements[0] = value_type(1);
    }

    /// @brief Sets all elements to random values in the [0, 1] interval.
    /// @details Only works if the Scalar type is not an int.
    void rand()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<value_type> dis(0.0, 1.0);
        m_elements[0] = dis(gen);
    }

    /// @brief Sets all elements to a specified value.
    void fill(value_type val)
    {
        m_elements[0] = val;
    }

    /// @brief Sets all elements using a lambda function.
    void fillf(std::function<value_type()> lambda)
    {
        m_elements[0] = lambda();
    }

    /// @brief Resets the vector to empty.
    void reset()
    {
        m_elements.clear();
    }

    /// @brief Returns the number of rows.
    size_t rows() const noexcept
    {
        return 1;
    }

    /// @brief Returns the number of columns.
    size_t cols() const noexcept
    {
        return 1;
    }

    /// @brief Printing method implementation, for printing the vector.
    void print(std::ostream& os) const
    {
        os << m_elements[0] << " ";
        os << "\n";
    }

private:
    storage_type m_elements;
};

}; // namespace lao

#endif // LAO_DENSE_MATRIX_H_