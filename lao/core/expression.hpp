#ifndef LAO_CORE_EXPRESSION_H_
#define LAO_CORE_EXPRESSION_H_

#include <cstddef>

namespace lao {

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

};

#endif // LAO_CORE_EXPRESSION_H_