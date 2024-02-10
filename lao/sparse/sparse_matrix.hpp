#ifndef LAO_SPARSE_MATRIX_H_
#define LAO_SPARSE_MATRIX_H_

#include <cstddef>
#include <lao/core/expression.hpp>

namespace lao {

/// @brief The class representing a sparse matrix.
/// @details This derives from the MatrixExpression class, and is a leaf in the
/// expression tree. This is a separate class to a dense matrix as it requires a
/// different storage mechanism and thus different access mechanism, due to
/// the sparsity of the matrix (being mostly 0s, and usually very large).
template <typename Scalar, size_t Rows, size_t Cols>
class SparseMatrix : public MatrixExpression<Scalar, SparseMatrix<Scalar, Rows, Cols>> {
public:
    /// @brief Default constructor which zero initializes the matrix.
    SparseMatrix()
        : m_elements(Rows * Cols, value_type(0))
    {
    }



    size_t rows() const { return Rows; }
    size_t cols() const { return Cols; }
};

}; // namespace lao

#endif // LAO_SPARSE_MATRIX_H_