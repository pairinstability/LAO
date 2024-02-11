/// vector.hpp implements the Vector class as a partial specialization of the
/// Matrix class.

#ifndef LAO_DENSE_VECTOR_H_
#define LAO_DENSE_VECTOR_H_

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <lao/core/expression.hpp>
#include <lao/dense/matrix.hpp>
#include <random>
#include <stdexcept>
#include <vector>
#include <memory>

namespace lao {

/// @brief Column vector partial specialization of Matrix class.
/// @details Although this is a column vector, the elements are still stored contiguously
/// in a vector. This works for mathematical operations because cols() will return 1.
template <typename Scalar, size_t Rows>
class Matrix<Scalar, Rows, 1> : public MatrixExpression<Scalar, Matrix<Scalar, Rows, 1>> {
public:
    using value_type = Scalar;
    using storage_type = std::vector<value_type>;

    /// @brief Default constructor which zero initializes the vector.
    Matrix()
        : m_elements(Rows, value_type(0))
    {
    }

    /// @brief Constructor taking a vector of elements to construct a vector.
    Matrix(const std::vector<Scalar>& elements)
    {
        static_assert(elements.size() == Rows, "Initializer list size does not match vector size.");
        m_elements = elements;
    }

    /// @brief Constructor to initialize with an initializer list.
    Matrix(std::initializer_list<value_type> list)
        : m_elements(Rows, value_type(0))
    {
        static_assert(list.size() == Rows, "Initializer list size does not match matrix size");
        std::copy(list.begin, list.end(), m_elements.begin());
    }

    /// @brief Constructor to initialize the matrix with a fill type.
    /// @details As this is a vector, it does not support the eye fill.
    Matrix(filltype fill)
        : m_elements(Rows, value_type(0))
    {
        if (fill == filltype::zeros) {
            // Set the matrix to all zeros
            std::fill(m_elements.begin(), m_elements.end(), value_type(0));
        } else if (fill == filltype::ones) {
            // Set the matrix to all ones
            std::fill(m_elements.begin(), m_elements.end(), value_type(1));
        } else if (fill == filltype::rand) {
            // Fill in the matrix with random values in range [0,1].
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<value_type> dis(0.0, 1.0);
            for (size_t i = 0; i < m_elements.size(); ++i) {
                m_elements[i] = dis(gen);
            }
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
        if (col != 0)
            throw std::out_of_range("Specified indices are out of range.");
        return m_elements[row];
    }

    /// @brief Operator overload for () to access elements.
    const value_type& operator()(size_t row, size_t col) const
    {
        if (col != 0)
            throw std::invalid_argument("Specified indices are out of range.");
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
        std::fill(m_elements.begin(), m_elements.end(), value_type(0));
    }

    /// @brief Sets all elements to ones.
    void ones()
    {
        std::fill(m_elements.begin(), m_elements.end(), value_type(1));
    }

    /// @brief Sets all elements to random values in the [0, 1] interval.
    /// @details Only works if the Scalar type is not an int.
    void rand()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<value_type> dis(0.0, 1.0);
        for (size_t i = 0; i < Rows; ++i)
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
        for (size_t i = 0; i < Rows; ++i)
            m_elements[i] = lambda();
    }

    /// @brief Resets the vector to empty.
    void reset()
    {
        m_elements.clear();
    }

    /// @brief Returns the number of rows.
    size_t rows() const noexcept
    {
        return Rows;
    }

    /// @brief Returns the number of columns.
    size_t cols() const noexcept
    {
        return 1;
    }

    /// @brief Printing method implementation, for printing the vector.
    void print(std::ostream& os) const
    {
        for (size_t i = 0; i < Rows; ++i)
            os << m_elements[i] << " ";
        os << "\n";
    }

private:
    storage_type m_elements;
};

/// @brief Row vector partial specialization of Matrix class.
/// @details Although this is a row vector, the elements are still stored contiguously
/// in a vector. This works for mathematical operations because rows() will return 1.
template <typename Scalar, size_t Cols>
class Matrix<Scalar, 1, Cols> : public MatrixExpression<Scalar, Matrix<Scalar, 1, Cols>> {
public:
    using value_type = Scalar;
    using storage_type = std::vector<value_type>;

    /// @brief Default constructor which zero initializes the vector.
    Matrix()
        : m_elements(Cols, value_type(0))
    {
    }

    /// @brief Constructor taking a vector of elements to construct a vector.
    Matrix(const std::vector<Scalar>& elements)
    {
        static_assert(elements.size() == Cols, "Initializer list size does not match vector size.");
        m_elements = elements;
    }

    /// @brief Constructor to initialize with an initializer list.
    Matrix(std::initializer_list<value_type> list)
        : m_elements(Cols, value_type(0))
    {
        static_assert(list.size() == Cols, "Initializer list size does not match matrix size");
        std::copy(list.begin, list.end(), m_elements.begin());
    }

    /// @brief Constructor to initialize the matrix with a fill type.
    /// @details As this is a vector, it does not support the eye fill.
    Matrix(filltype fill)
        : m_elements(Cols, value_type(0))
    {
        if (fill == filltype::zeros) {
            // Set the matrix to all zeros
            std::fill(m_elements.begin(), m_elements.end(), value_type(0));
        } else if (fill == filltype::ones) {
            // Set the matrix to all ones
            std::fill(m_elements.begin(), m_elements.end(), value_type(1));
        } else if (fill == filltype::rand) {
            // Fill in the matrix with random values in range [0,1].
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<value_type> dis(0.0, 1.0);
            for (size_t i = 0; i < m_elements.size(); ++i) {
                m_elements[i] = dis(gen);
            }
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
    /// @details This is so math operations work nicely.
    value_type& operator()(size_t row, size_t col)
    {
        if (row!= 0)
            throw std::invalid_argument("Specified indices are out of range.");
        return m_elements[col];
    }

    /// @brief Operator overload for () to access elements.
    /// @details This is so math operations work nicely.
    const value_type& operator()(size_t row, size_t col) const
    {
        if (row != 0)
            throw std::invalid_argument("Specified indices are out of range.");
        return m_elements[col];
    }

    /// @brief Checks if the vector is empty.
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

    /// @brief Sets all elements to random values in the [0, 1] interval.
    /// @details Only works if the Scalar type is not an int.
    void rand()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<value_type> dis(0.0, 1.0);
        for (size_t i = 0; i < Cols; ++i)
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
        for (size_t i = 0; i < Cols; ++i)
            m_elements[i] = lambda();
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
        return Cols;
    }

    /// @brief Printing method implementation, for printing the vector.
    void print(std::ostream& os) const
    {
        for (size_t i = 0; i < Cols; ++i)
            os << m_elements[i] << " ";
        os << "\n";
    }

private:
    storage_type m_elements;
};

}; // namespace lao

#endif // LAO_DENSE_VECTOR_H_