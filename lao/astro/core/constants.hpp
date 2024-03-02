#ifndef LAO_ASTRO_CORE_CONSTANTS_H_
#define LAO_ASTRO_CORE_CONSTANTS_H_

#include <array>
#include <iostream>

namespace lao {
namespace astro {

#define MAKE_CONSTANT(symbol, value) \
    template <typename T>            \
    static constexpr T symbol = value;

    // pi, dimensionless
    MAKE_CONSTANT(PI, 3.141592653589793238462643383279502884);
    // astronomical unit [m]
    MAKE_CONSTANT(AU, 1.4959787070691e11);
    // speed of light [m/s]
    MAKE_CONSTANT(C, 299792458);
    // gravity [m/s^2]
    MAKE_CONSTANT(G, 9.80665);

    // standard gravitional parameters [m^3 / s^-2]
    MAKE_CONSTANT(MU_SUN, 1.327124400189e20);
    MAKE_CONSTANT(MU_MERCURY, 2.20329e13);
    MAKE_CONSTANT(MU_VENUS, 3.248599e14);
    MAKE_CONSTANT(MU_EARTH, 3.9860044188e14);
    MAKE_CONSTANT(MU_MOON, 4.90486959e12);
    MAKE_CONSTANT(MU_MARS, 4.2828372e13);
    MAKE_CONSTANT(MU_JUPITER, 1.266865349e17);
    MAKE_CONSTANT(MU_SATURN, 3.79311879e16);
    MAKE_CONSTANT(MU_URANUS, 5.7939399e15);
    MAKE_CONSTANT(MU_NEPTUNE, 6.8365299e15);
    MAKE_CONSTANT(MU_PLUTO, 8.719e11);

    // earth radius [m]
    MAKE_CONSTANT(EARTH_RADIUS, 6.3781366e6)

    // degrees to radians
    MAKE_CONSTANT(DEG2RAD, PI<T> / 180.0);
    // radians to degrees
    MAKE_CONSTANT(RAD2DEG, 180.0 / PI<T>);
    // days to seconds
    MAKE_CONSTANT(DAY2SEC, 86400.0);
    // seconds to day
    MAKE_CONSTANT(SEC2DAY, 1. / 86400.0);
    // day to year
    MAKE_CONSTANT(DAY2YEAR, 1. / 365.25);
    // AU to m
    MAKE_CONSTANT(AU2M, 149597870691.0);

    // solver tolerance
    MAKE_CONSTANT(SOLVER_TOLERANCE, 1e-16);

}; // namespace astro
}; // namespace lao

#endif // LAO_ASTRO_CORE_CONSTANTS_H_