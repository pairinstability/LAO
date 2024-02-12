#ifndef LAO_SPARSE_MATRIX_H_
#define LAO_SPARSE_MATRIX_H_

#include <cstddef>
#include <fstream>
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
        csr.m_rowvec.resize(R + 1, value_type(0));
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
                    csr.insert(row, col, val);
                ++col;
            }
            ++row;
        }
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

    CSRStorage csr;
};

}; // namespace lao

#endif // LAO_SPARSE_MATRIX_H_