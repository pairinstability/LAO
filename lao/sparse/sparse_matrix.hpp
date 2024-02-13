#ifndef LAO_SPARSE_MATRIX_H_
#define LAO_SPARSE_MATRIX_H_

#include <cstddef>
#include <fstream>
#include <functional>
#include <lao/core/expression.hpp>
#include <lao/core/forward.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lao {

/// @brief The class representing a sparse matrix.
/// @details A matrix is represented by a number of rows and columns,
/// noted as `row x column`.
/// It is templated with a Scalar parameter, and a Row and Column.
/// The storage mechanism is using std::vectors with compressed sparse row format.
template <typename S, size_t R, size_t C>
class SparseMatrix : public MatrixExpression<SparseMatrix<S, R, C>, S, R, C> {
public:
    using value_type = S;
    using storage_type_v = std::vector<value_type>;
    using storage_type_row = std::vector<size_t>;
    using storage_type_col = std::vector<size_t>;

    SparseMatrix()
    {
        m_csr.m_rowvec.resize(R + 1, value_type(0));
    }

    SparseMatrix(const std::string& filename)
        : SparseMatrix()
    {
        std::ifstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Failed to open file");

        std::string line;
        size_t row = 0;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string cell;
            size_t col = 0;
            while (std::getline(iss, cell, ',')) {
                value_type val = std::stod(cell);
                if (val != 0)
                    m_csr.insert(row, col, val);
                ++col;
            }
            ++row;
        }
    }

    /// @brief Copy constructor.
    SparseMatrix(const SparseMatrix& other)
        : m_csr(other.m_csr)
    {
    }

    SparseMatrix& operator=(const SparseMatrix& other)
    {
        if (&other != this)
            m_csr = other.m_csr;
        return *this;
    }

    /// @brief Operator for converting MatrixExpression <-> SparseMatrix
    template <typename E>
    SparseMatrix(const MatrixExpression<E, S, R, C>& expr)
        : SparseMatrix()
    {
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < C; ++j) {
                value_type val = static_cast<value_type>(expr(i, j));
                if (val != 0)
                    m_csr.insert(i, j, val);
            }
        }
    }

    /// @brief operator overload for () to access elements.
    value_type& operator()(size_t row, size_t col)
    {
        if (row > R || col > C)
            throw std::out_of_range("Specified indices are out of range.");
        return m_csr.get(row, col);
    }

    /// @brief operator overload for () to access elements.
    const value_type& operator()(size_t row, size_t col) const
    {
        if (row > R || col > C)
            throw std::out_of_range("Specified indices are out of range.");
        return m_csr.get(row, col);
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
        return m_csr.m_values.empty();
    }

    /// @brief Sets all elements to zero.
    /// @details Given this is a sparse matrix, we only store non-zero values so this
    /// clears all of the storage vectors used in CSR.
    void zeros()
    {
        m_csr.m_values.clear();
        m_csr.m_colvec.clear();
        m_csr.m_rowvec.clear();
    }

    /// @brief Sets all elements to the identity matrix.
    void eye()
    {
        if constexpr (R == C) {
            zeros();
            for (size_t i = 0; i < std::min(R, C); ++i)
                m_csr.insert(i, i, 1);
        } else {
            throw std::logic_error("Identity matrix is only defined for square matrices.");
        }
    }

    /// @brief Sets all elements using a lambda function.
    /// @details This is applying to every element, not just the non-zero elements.
    void fillf(std::function<value_type()> lambda)
    {
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < C; ++j) {
                value_type val = m_csr.get(i, j);
                if (val != 0) {
                    m_csr.insert(i, j, lambda());
                }
            }
        }
    }

    /// @brief Sets only the non-zero elements using a lambda function.
    void fillfnz(std::function<value_type()> lambda)
    {
        for (size_t i = 0; i < m_csr.rows(); ++i) {
            for (auto it = m_csr.row_begin(i); it != m_csr.row_end(i); ++it) {
                size_t j = it->column;
                value_type val = it->value;
                m_csr.insert(i, j, lambda());
            }
        }
    }

    /// @brief Resets the matrix to empty.
    /// @details Due to how the sparse matrix storage only stores non-zero values,
    /// this is functionally equivalent to just zero-ing the matrix.
    void reset()
    {
        zeros();
    }

private:
    struct CSRStorage {
        // Compressed sparse row format uses three vectors to store information about
        // a sparse matrix.
        // Let NNZ be the number of non-zero values.
        // contains the non-zero values. size NNZ.
        storage_type_v m_values;
        // contains row indices. size m + 1. The last element is NNZ.
        storage_type_row m_rowvec;
        // contains the column indices. size NNZ.
        storage_type_col m_colvec;

        void insert(size_t row, size_t col, value_type val)
        {
            m_values.push_back(val);
            m_colvec.push_back(col);
            for (size_t i = row + 1; i < m_rowvec.size(); ++i) {
                ++m_rowvec[i];
            }
        }

        value_type get(size_t row, size_t col) const
        {
            size_t start = m_rowvec[row];
            size_t end = m_rowvec[row + 1];
            for (size_t i = start; i < end; ++i) {
                if (m_colvec[i] == col) {
                    return m_values[i];
                }
            }
            return value_type(0);
        }
    };

    CSRStorage m_csr;
};

}; // namespace lao

#endif // LAO_SPARSE_MATRIX_H_