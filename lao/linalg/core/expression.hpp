#ifndef LAO_LINALG_CORE_EXPRESSION_H_
#define LAO_LINALG_CORE_EXPRESSION_H_

#include <cstddef>

namespace lao {
namespace linalg {

    template <typename Derived, typename S, size_t R, size_t C>
    class MatrixExpression {
    public:
        S operator()(size_t row, size_t col) const
        {
            return static_cast<const Derived*>(this)->operator()(row, col);
        }

        size_t rows() const noexcept
        {
            return R;
        }

        size_t cols() const noexcept
        {
            return C;
        }
    };

}; // namespace linalg
}; // namespace lao

#endif // LAO_LINALG_CORE_EXPRESSION_H_