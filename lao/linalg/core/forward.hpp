#ifndef LAO_LINALG_CORE_FORWARD_H_
#define LAO_LINALG_CORE_FORWARD_H_

#include <cstddef>
#include <vector>

namespace lao {
namespace linalg {

    template <typename S, size_t R, size_t C, typename B = std::vector<S>>
    class Matrix;

    template <typename S, size_t R, size_t C>
    class SparseMatrix;

    template <typename Derived, typename S, size_t R, size_t C>
    class MatrixExpression;

}; // namespace linalg
}; // namespace lao

#endif // LAO_LINALG_CORE_FORWARD_H_