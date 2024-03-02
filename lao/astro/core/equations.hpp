#ifndef LAO_ASTRO_EQUATIONS_H_
#define LAO_ASTRO_EQUATIONS_H_

#include <cmath>

namespace lao {
namespace astro {

    /// @brief equation representing the conversion from eccentric anomaly and eccentricity to mean anomaly
    /// @details https://en.wikipedia.org/wiki/Mean_anomaly
    /// @param E eccentric anomaly
    /// @param e eccentricity
    /// @param M mean anomaly
    /// @returns
    inline double meanAnomaly(double E, double e, double M)
    {
        return E - e * std::sin(E) - M;
    }

    /// @brief equation representing the derivative of the conversion from eccentric anomaly and eccentricity to mean anomaly
    /// @details https://en.wikipedia.org/wiki/Mean_anomaly
    /// @param E eccentric anomaly
    /// @param e eccentricity
    /// @returns
    inline double meanAnomalyDerivative(double E, double e)
    {
        return 1.0 - e * std::cos(E);
    }

}; // namespace astro
}; // namespace lao

#endif // KEPELR_EQUATIONS_H_