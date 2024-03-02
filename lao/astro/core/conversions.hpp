#ifndef LAO_ASTRO_CONVERSIONS_H_
#define LAO_ASTRO_CONVERSIONS_H_

#include <cmath>
#include <functional>
#include <lao/astro/core/constants.hpp>
#include <lao/astro/core/equations.hpp>
#include <lao/astro/math/newton_raphson.hpp>
#include <lao/astro/math/transforms.hpp>
#include <lao/linalg/dense/matrix.hpp>

namespace lao {
namespace astro {

    /// @brief converts the mean anomaly M to eccentric anomaly E.
    /// @details uses the newton-raphson method using the mean anomaly equation to solve for E. see
    /// https://ssd.jpl.nasa.gov/planets/approx_pos.html.
    /// @param M mean anomaly.
    /// @param e eccentricity.
    /// @returns eccentric anomaly.
    template <typename T>
    inline T meanAnomalyToEccentricAnomaly(const T& M, const T& e)
    {
        T E = M + e * sin(M);

        auto F = [&](T E) { return meanAnomaly(E, e, M); };
        auto dF = [&](T E) { return meanAnomalyDerivative(E, e); };

        return newtonRaphson(E, F, dF, 100, SOLVER_TOLERANCE<T>);
    };

    /// @brief converts keplerian elements to cartesian coordinate state vector.
    /// @details see https://ssd.jpl.nasa.gov/planets/approx_pos.html.
    /// @param elements keplerian elements as a 6D vector [a,e,i,Omega,omega,E].
    /// @param mu_central_body the standard gravitational parameter of the attracting body.
    /// @returns the cartesian state vector as a 6D vector.
    template <typename T>
    inline linalg::RowVector<double, 6> keplerianToCartesian(const linalg::RowVector<double, 6>& elements, const T& mu_central_body)
    {
        T a = elements(1, 1);
        T e = elements(1, 2);
        T i = elements(1, 3);
        T Omega = elements(1, 4);
        T omega = elements(1, 5);
        T E = elements(1, 6);

        linalg::RowVector<double, 3> angles;

        angles(1, 1) = -omega;
        angles(1, 2) = -i;
        angles(1, 3) = -Omega;

        // computing the body's perifocal coordinates in its orbital plane
        // https://en.wikipedia.org/wiki/Perifocal_coordinate_system
        // https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html#step-1transform-to-perifocal-frame
        linalg::Matrix<T, 3, 3> R = rotationFromEuler("ZXZ", angles);

        // position and velocity vectors in the perifocal frame
        T xprime = a * (std::cos(E) - e);
        T yprime = a * (std::sqrt(1 - std::pow(e, 2))) * std::sin(E);
        T zprime = 0;
        T vxprime = -std::sqrt(mu_central_body / a) * std::sin(E);
        T vyprime = std::sqrt(mu_central_body / a) * std::sqrt(1 - e * e) * std::cos(E);
        T vzprime = 0;

        linalg::Matrix<T, 1, 3> rprime;
        linalg::Matrix<T, 1, 3> vprime;

        rprime(1, 1) = xprime;
        rprime(1, 2) = yprime;
        rprime(1, 3) = zprime;
        vprime(1, 1) = vxprime;
        vprime(1, 2) = vyprime;
        vprime(1, 3) = vzprime;

        // 1x3 = 1x3 * 3x3
        linalg::Matrix<T, 1, 3> rvec = rprime * R;
        linalg::Matrix<T, 1, 3> vvec = vprime * R;

        linalg::Matrix<T, 1, 6> result = concat(rvec, vvec);

        return result;
    }

}; // namespace astro
}; // namespace lao

#endif // LAO_ASTRO_CONVERSIONS_H_