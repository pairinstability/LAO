/// this contains the base class and associated features for defining a planetary body.
/// a body can be described by:
///
/// * standard gravitational parameter [m^3/s^2]
/// * standard gravitational parameter of its central/parent body (assuming a 2-body system) [m^3/s^2]
/// * body radius [m]
/// * mass of the body [kg]
///
/// at any given moment in time, this body can also be described by:
/// * 3D cartesian position vector [m, m, m]
/// * 3D cartesian velocity vector [m/s, m/s, m/s]
///
/// which can be used to describe the orbital elements rleative to the parent body:
/// * eccentricity (e)
/// * semi-major axis (a)
/// * inclination (i)
/// * longitude of the ascending node (Omega)
/// * argument of peripapsis (omega)
/// * mean anomaly (M)

#ifndef LAO_ASTRO_BODY_BASE_H_
#define LAO_ASTRO_BODY_BASE_H_

#include <lao/astro/core/constants.hpp>
#include <lao/astro/date/epoch.hpp>
#include <lao/linalg/core/forward.hpp>
#include <lao/linalg/dense/matrix.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lao {
namespace astro {

    /// @brief parent class that represents a body.
    /// @details designed to be inherited such that the child implements the eph method
    /// which returns the state vector of the body using a child-specific method.
    class Base {
    protected:
        double m_mu_body;
        double m_mu_central_body;
        double m_radius;
        std::string m_name;

    public:
        /// @brief Base constructor.
        /// @param mu_body standard gravitational parameter of the body.
        /// @param mu_central_body standard gravitational paramater of the attracting/parent body.
        /// @param name name of the body.
        Base(double mu_body, double mu_central_body, double radius, const std::string& name)
            : m_mu_body(mu_body)
            , m_mu_central_body(mu_central_body)
            , m_radius(radius)
            , m_name(name)
        {
            if (radius <= 0)
                throw std::invalid_argument("Radius must be greater than zero");
            if (mu_body <= 0)
                throw std::invalid_argument("Body standard gravitational body must be greater than zero");
            if (mu_central_body <= 0)
                throw std::invalid_argument("Central body standard gravitational body must be greater than zero");
        };

        virtual ~Base() {};

        /// @brief returns the cartesian coordinate form of position and velocity vectors given epoch.
        /// @details this is a pure virtual method which is implemented in inheriting classes.
        /// @param epoch_date the epoch from which to find the position and velocity of the body.
        /// @returns a 6D array representing the position and velocity vectors.
        virtual linalg::RowVector<double, 6> eph(const Epoch& epoch_date) const = 0;

        virtual std::string ostreamExtra() const
        {
            return std::string();
        }

        friend std::ostream& operator<<(std::ostream& os, const Base& body)
        {
            os << "{\n";
            os << "\"body\": \"" << body.m_name << "\",\n";
            os << "\"gravitational_parameter_m3_per_s2\": " << body.m_mu_body << ",\n";
            os << "\"parent_gravitational_parameter_m3_per_s2\": " << body.m_mu_central_body << ",\n";
            os << "\"body_radius_m\": " << body.m_radius << ",\n";
            os << body.ostreamExtra();
            os << "}\n";
            return os;
        };
    };

}; // namespace astro 
}; // namespace lao

#endif // LAO_ASTRO_BODY_BASE_H_