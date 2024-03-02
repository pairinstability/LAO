/// JPL low precision uses a lower accuracy formula for planetary positions: https://ssd.jpl.nasa.gov/planets/approx_pos.html

#ifndef LAO_ASTRO_BODY_JPL_LOW_PRECISION_H_
#define LAO_ASTRO_BODY_JPL_LOW_PRECISION_H_

#include <lao/astro/body/base.hpp>
#include <lao/astro/core/constants.hpp>
#include <lao/astro/core/conversions.hpp>
#include <lao/linalg/dense/matrix.hpp>
#include <lao/astro/date/epoch.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace lao {
namespace astro {

    /// @brief JPL low precision class
    /// @details inherits the Base class representing an abstract body. this class implements
    /// the JPL low precision ephemerides. it uses a lookup hashtable and given a name of the body
    /// in the constructor to find the parameters and orbital parameters to calculate the ephemerides.
    class JPLLP : public Base {
    private:
        // table of bodies and some basic parameters, taken from numerous sources.
        //
        // +----------+-----------------+---------------------+------------+
        // |          |   mu [m^3/s^2]  | mu parent [m^3/s^2] | radius [m] |
        // | mercury  | 2.20329e13      | 1.327124400189e20   | 2439500    |
        // | venus    | 3.248599e14     | 1.327124400189e20   | 6052000    |
        // | E-M bary | 3.9860044188e14 | 1.327124400189e20   | 6378100    |
        // | mars     | 4.2828372e13    | 1.327124400189e20   | 3396000    |
        // | jupiter  | 1.266865349e17  | 1.327124400189e20   | 71492000   |
        // | saturn   | 3.79311879e18   | 1.327124400189e20   | 60268000   |
        // | uranus   | 5.7939399e15    | 1.327124400189e20   | 25559000   |
        // | neptune  | 6.8365299e15    | 1.327124400189e20   | 24764000   |
        // +----------+-----------------+---------------------+------------+
        //
        // taken from https://ssd.jpl.nasa.gov/planets/approx_pos.html
        // keplerian elements and their rates, with respect to the mean ecliptic and equinox of J2000,
        // valid for the time-interval 1800 AD - 2050 AD.
        //
        // +------------+-------------+------------+-------------+--------------+-----------------+-----------------+
        // |            |    a [au]   |   e [rad]  |   I [deg]   |    L [deg]   | long peri [deg] | long node [deg] |
        // | mercury    | 0.38709927  | 0.20563593 | 7.00497902  | 252.25032350 | 77.45779628     | 48.33076593     |
        // | venus      | 0.72333566  | 0.00677672 | 3.39467605  | 181.97909950 | 131.60246718    | 76.67984255     |
        // | E-M bary   | 1.00000261  | 0.01671123 | -0.00001531 | 100.46457166 | 102.93768193    | 0.0             |
        // | mars       | 1.52371034  | 0.09339410 | 1.84969142  | -4.55343205  | -23.94362959    | 49.55953891     |
        // | jupiter    | 5.20288700  | 0.04838624 | 1.30439695  | 34.39644051  | 14.72847983     | 100.47390909    |
        // | saturn     | 9.53667594  | 0.05386179 | 2.48599187  | 49.95424423  | 92.59887831     | 113.66242448    |
        // | uranus     | 19.18916464 | 0.04725744 | 0.77263783  | 313.23810451 | 170.95427630    | 74.01692503     |
        // | neptune    | 30.06992276 | 0.00859048 | 1.77004347  | -55.12002969 | 44.96476227     | 131.78422574    |
        // +------------+-------------+------------+-------------+--------------+-----------------+-----------------+
        //
        // +------------+--------------+---------------+---------------+-----------------+-----------------------+-----------------------+
        // |            | adot [au/Cy] | edot [rad/Cy] | Idot [deg/Cy] | Ldot [deg/Cy]   | long peridot [deg/Cy] | long nodedot [deg/Cy] |
        // | mercury    | 0.00000037   | 0.00001906    | -0.00594749   | 149472.67411175 | 0.16047689            | -0.12534081           |
        // | venus      | 0.00000390   | -0.00004107   | -0.00078890   | 58517.81538729  | 0.00268329            | -0.27769418           |
        // | E-M bary   | 0.00000562   | -0.00004392   | -0.01294668   | 35999.37244981  | 0.32327364            | 0.0                   |
        // | mars       | 0.00001847   | 0.00007882    | -0.00813131   | 19140.30268499  | 0.44441088            | -0.29257343           |
        // | jupiter    | -0.00011607  | -0.00013253   | -0.00183714   | 3034.74612775   | 0.21252668            | 0.20469106            |
        // | saturn     | -0.00125060  | -0.00050991   | 0.00193609    | 1222.49362201   | -0.41897216           | -0.28867794           |
        // | uranus     | -0.00196176  | -0.00004397   | -0.00242939   | 428.48202785    | 0.40805281            | 0.04240589            |
        // | neptune    | 0.00026291   | 0.00005105    | 0.00035372    | 218.45945325    | -0.32241464           | -0.00508664           |
        // +------------+--------------+---------------+---------------+-----------------+-----------------------+-----------------------+
        //
        // although we're currently only using solar system large bodies and linear ordering means we could just use
        // array indexing, this is expandable for any number of bodies and a hashmap will allow for fast lookup
        //
        // the table is ordered as such: [body parameters x3], [keplerian elements x6], [keplerian element rates x6]
        // the table below contains the bodies and some basic parameters, and then their keplerian elements.

        struct Body {
            double mu;
            double mu_central;
            double radius;
            double a;
            double e;
            double I;
            double L;
            double long_peri;
            double long_node;
            double adot;
            double edot;
            double Idot;
            double Ldot;
            double long_peridot;
            double long_nodedot;
        };

        // TODO: perf this against a sorted vector since with a smaller number of elements it may be better
        // TODO: use the 3000 BC - 3000 AD table rather than this one
        std::unordered_map<std::string, Body> m_body_lookup_1800_2050 = {
            { "Mercury", { MU_MERCURY<double>, MU_SUN<double>, 2439500, 0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593, 0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081 } },
            { "Venus", { MU_VENUS<double>, MU_SUN<double>, 6052000, 0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255, 0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418 } },
            { "EM bary", { MU_EARTH<double>, MU_SUN<double>, 6378100, 1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0, 0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0 } },
            { "Mars", { MU_MARS<double>, MU_SUN<double>, 3396000, 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891, 0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343 } },
            { "Jupiter", { MU_JUPITER<double>, MU_SUN<double>, 71492000, 5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909, -0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106 } },
            { "Saturn", { MU_SATURN<double>, MU_SUN<double>, 60268000, 9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448, -0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794 } },
            { "Uranus", { MU_URANUS<double>, MU_SUN<double>, 25559000, 19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503, -0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589 } },
            { "Neptune", { MU_NEPTUNE<double>, MU_SUN<double>, 24764000, 30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574, 0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664 } }
        };

        // JPL LP elements position [x,y,z]
        linalg::RowVector<double, 6> m_jpl_elements;
        // JPL LP elements velocity [xdot, ydot, zdot]
        linalg::RowVector<double, 6> m_jpl_elements_dot;

    public:
        /// @brief JPL low precision constructor.
        /// @param body_name name of the body
        /// TODO: better initial arguments
        JPLLP(const std::string& body_name)
            : Base(1.0, 1.0, 1.0, body_name)
        {
            // perform a lookup and set the corresponding values given the body name
            if (auto it = m_body_lookup_1800_2050.find(body_name); it != m_body_lookup_1800_2050.end()) {
                m_mu_body = it->second.mu;
                m_mu_central_body = it->second.mu_central;
                m_radius = it->second.radius;

                m_jpl_elements(1,1) = it->second.a;
                m_jpl_elements(1,2) = it->second.e;
                m_jpl_elements(1,3) = it->second.I;
                m_jpl_elements(1,4) = it->second.L;
                m_jpl_elements(1,5) = it->second.long_peri;
                m_jpl_elements(1,6) = it->second.long_node;

                m_jpl_elements_dot(1,1) = it->second.adot;
                m_jpl_elements_dot(1,2) = it->second.edot;
                m_jpl_elements_dot(1,3) = it->second.Idot;
                m_jpl_elements_dot(1,4) = it->second.Ldot;
                m_jpl_elements_dot(1,5) = it->second.long_peridot;
                m_jpl_elements_dot(1,6) = it->second.long_nodedot;

                //                m_jpl_elements = { it->second.a, it->second.e, it->second.I, it->second.L, it->second.long_peri, it->second.long_node };
                //                m_jpl_elements_dot = { it->second.adot, it->second.edot, it->second.Idot, it->second.Ldot, it->second.long_peridot, it->second.long_nodedot };
            } else {
                throw std::invalid_argument("Unknown body name");
            }
        };

        /// @brief returns the cartesian coordinate form of position and velocity vectors given epoch.
        /// @details this uses the JPL low precision ephemerides.
        /// @param epoch_date the epoch from which to find the position and velocity of the body.
        /// @returns a 6D array representing the position and velocity vectors.
        linalg::RowVector<double, 6> eph(const Epoch& epoch_date) const
        {
            // converting from elements to pos and vel using https://ssd.jpl.nasa.gov/planets/approx_pos.html
            /// TODO: se the 3000 BC - 3000 AD table rather than this one
            if (epoch_date.MJD2000() <= -73048.0 || epoch_date.MJD2000() >= 18263.0)
                throw std::invalid_argument("epoch date must be in range [1800, 2050]");

            // using equation from https://ssd.jpl.nasa.gov/planets/approx_pos.html to obtain coordinates from ephemeris
            // 1. compute the value of each of that planet's six elements, a = a0 + adot * T where T is the number of centuries past J2000.0
            // 2. convert units to SI
            // 3. compute argument of periapsis and mean anomaly
            // 4. solve equation
            // 5. convert keplerian elements (a,e,i,W,w,E) to cartesian

            // elements contains:
            // * semi-major axis, a [m]
            // * eccentricity, e
            // * inclination, i [rad]
            // * longitude of the ascending node, Omega [rad]
            // * argument of periapsis, omega [rad]
            // * mean anomaly, M
            linalg::RowVector<double, 6> elements;

            // 1.
            // number of centuries past J2000.0
            double dt = (epoch_date.MJD2000() - 2451545.0) / 36525.0;
            // a = a0 + adot * T and converting to SI
            // a [AU -> m]
            // https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes
            elements(1,1) = (m_jpl_elements(1,1) + m_jpl_elements_dot(1,1) * dt) * AU2M<double>;
            // e
            // https://en.wikipedia.org/wiki/Orbital_eccentricity
            elements(1,2) = (m_jpl_elements(1,2) + m_jpl_elements_dot(1,2) * dt);
            // i [deg -> rad]
            // https://en.wikipedia.org/wiki/Orbital_inclination
            elements(1,3) = (m_jpl_elements(1,3) + m_jpl_elements_dot(1,3) * dt) * DEG2RAD<double>;
            // Omega [deg -> rad]
            // https://en.wikipedia.org/wiki/Longitude_of_the_ascending_node
            elements(1,4) = (m_jpl_elements(1,4) + m_jpl_elements_dot(1,4) * dt) * DEG2RAD<double>;
            // omega [deg -> rad]; omega = omegbar - Omega
            // https://en.wikipedia.org/wiki/Argument_of_periapsis
            elements(1,5) = (m_jpl_elements(1,5) + m_jpl_elements_dot(1,5) * dt) * DEG2RAD<double> - elements(1,4);
            // 4.
            // M
            elements(1,6) = meanAnomalyToEccentricAnomaly<double>(elements(1,6), elements(1,2));
            // 5.
            return keplerianToCartesian<double>(elements, m_mu_central_body);
        };

        // TODO: consistent accuracy
        std::string ostreamExtra() const
        {
            std::ostringstream s;
            s << "\"JPL_low_precision\": {\n";
            s << "  \"semi_major_axis_au\": " << m_jpl_elements(1,1) << ",\n";
            s << "  \"eccentricity\": " << m_jpl_elements(1,2) << ",\n";
            s << "  \"inclination_deg\": " << m_jpl_elements(1,3) << ",\n";
            s << "  \"mean_longitude_deg\": " << m_jpl_elements(1,4) << ",\n";
            s << "  \"longitude_of_perihelion_deg\": " << m_jpl_elements(1,5) << ",\n";
            s << "  \"longitude_of_ascending_node_deg\": " << m_jpl_elements(1,6) << ",\n";
            s << "  \"semi_major_axis_rate_of_change_au_per_Cy\": " << m_jpl_elements_dot(1,1) << ",\n";
            s << "  \"eccentricity_rate_of_change_rad_per_Cy\": " << m_jpl_elements_dot(1,2) << ",\n";
            s << "  \"inclination_rate_of_change_deg_per_Cy\": " << m_jpl_elements_dot(1,3) << ",\n";
            s << "  \"mean_longitude_rate_of_change_deg_per_Cy\": " << m_jpl_elements_dot(1,4) << ",\n";
            s << "  \"longitude_of_perihelion_rate_of_change_deg_per_Cy\": " << m_jpl_elements_dot(1,5) << ",\n";
            s << "  \"longitude_of_ascending_node_rate_of_change_deg_per_Cy\": " << m_jpl_elements_dot(1,6) << "\n";
            s << "  }\n";
            return s.str();
        }
    };

}; // namespace astro
}; // namespace lao

#endif // LAO_ASTRO_BODY_JPL_LOW_PRECISION_H_