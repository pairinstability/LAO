#ifndef LAO_LINALG_MATH_CONSTRAINTS_H_
#define LAO_LINALG_MATH_CONSTRAINTS_H_

#include <cstddef>
#include <functional>

namespace lao {
namespace linalg {

    ///// @brief Concept for matrix invertability.
    // template <typename S, size_t R, size_t C>
    // concept EnforceInvertability = requires
    //{
    // };

    /// @brief Concept for enforcing same scalar type.
    template <typename S1, typename S2>
    concept EnforceSameType = requires
    {
        requires std::same_as<S1, S2>;
    };

    /// @brief Concept for enforcing matrix addition/subtraction constraints.
    template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2>
    concept EnforceSameShape = requires
    {
        // Check if the scalar types match
        requires std::same_as<S1, S2>;
        // Check if the dimensions match
        requires R1 == R2&& C1 == C2;
    };

    /// @brief Concept for enforcing matrix multiplication constraints.
    template <typename S1, typename S2, size_t R1, size_t R2, size_t C1, size_t C2>
    concept EnforceMatMulReqs = requires
    {
        // Check if the scalar types match
        requires std::same_as<S1, S2>;
        // Check if the dimensions match
        requires C1 == R2;
    };

    /// @brief Concept for enforcing square matrices.
    template <typename S, size_t R, size_t C>
    concept EnforceSquareMatrix = requires
    {
        requires R == C;
    };

}; // namespace linalg
}; // namespace lao

#endif // LAO_LINALG_MATH_CONSTRAINTS_H_