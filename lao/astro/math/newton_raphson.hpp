/// solving the newton-raphson method, or newton's method for finding approximations for roots
/// https://en.wikipedia.org/wiki/Newton%27s_method

#ifndef LAO_ASTRO_MATH_NEWTON_RAPHSON_H_
#define LAO_ASTRO_MATH_NEWTON_RAPHSON_H_

#include <algorithm>
#include <cmath>

namespace lao {
namespace astro {

    /// @brief generically solving the newton-raphson method for finding approximations for roots.
    /// @param x starting point.
    /// @param F function.
    /// @param dF derivative of said function.
    /// @param max_iterations maximum number of loop iterations.
    /// @param accuracy accuracy.
    /// @returns solution.
    template <class start, class func, class deriv>
    inline double newtonRaphson(start& x, func F, deriv dF, double max_iterations, const double& accuracy)
    {
        start term;
        do {
            term = F(x) / dF(x);
            x -= term;
        } while ((std::fabs(term / std::max(std::fabs(x), 1.0)) > accuracy) && (--max_iterations));
        return max_iterations;
    };

}; // namespace astro
}; // namespace lao

#endif // LAO_ASTRO_MATH_NEWTON_RAPHSON_H_