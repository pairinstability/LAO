#ifndef LAO_ASTRO_MATH_LINEAR_TRANSFORMS_H_
#define LAO_ASTRO_MATH_LINEAR_TRANSFORMS_H_

#include <cmath>
#include <lao/linalg/dense/matrix.hpp>

namespace lao {
namespace astro {

    /// @brief builds a rotation matrix using euler angles.
    /// @details see https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations.
    /// @param axis axis from which the 3D rotation is around
    /// @param angles euler angles [rad,rad,rad]
    template <typename T>
    linalg::Matrix<T, 3, 3> rotationFromEuler(const std::string& axis, const linalg::RowVector<T, 3>& angles)
    {
        T phi = angles(1, 1);
        T theta = angles(1, 2);
        T psi = angles(1, 3);

        linalg::Matrix<T, 3, 3> elements;

        float c1 = std::cos(phi);
        float c2 = std::cos(theta);
        float c3 = std::cos(psi);
        float s1 = std::sin(phi);
        float s2 = std::sin(theta);
        float s3 = std::sin(psi);

        if (axis == "ZXZ") {
            elements(1, 1) = c2;
            elements(1, 2) = -c3 * s2;
            elements(1, 3) = s2 * s3;
            elements(2, 1) = c1 * s2;
            elements(2, 2) = c1 * c2 * c3 - s1 * s3;
            elements(2, 3) = -c3 * s1 - c1 * c2 * s3;
            elements(3, 1) = s1 * s2;
            elements(3, 2) = c1 * s3 + c2 * c3 * s1;
            elements(3, 3) = c1 * c3 - c2 * s1 * s3;
        } else if (axis == "XYX") {
            elements(1, 1) = c2;
            elements(1, 2) = s2 * s3;
            elements(1, 3) = c3 * s2;
            elements(2, 1) = s1 * s2;
            elements(2, 2) = c1 * c3 - c2 * s1 * s3;
            elements(2, 3) = -c1 * s3 - c2 * c3 * s1;
            elements(3, 1) = -c1 * s2;
            elements(3, 2) = c3 * s1 + c1 * c2 * s3;
            elements(3, 3) = c1 * c2 * c3 - s1 * s3;
        } else if (axis == "YXY") {
            elements(1, 1) = c1 * c3 - c2 * s1 * s3;
            elements(1, 2) = s1 * s2;
            elements(1, 3) = c1 * s3 + c2 * c3 * s1;
            elements(2, 1) = c2 * s3;
            elements(2, 2) = c2;
            elements(2, 3) = -c3 * s2;
            elements(3, 1) = -c3 * s1 - c1 * c2 * s3;
            elements(3, 2) = c1 * s2;
            elements(3, 3) = c1 * c2 * c3 - s1 * s3;
        } else if (axis == "YZY") {
            elements(1, 1) = c1 * c2 * c3 - s1 * s3;
            elements(1, 2) = -c1 * s2;
            elements(1, 3) = c3 * s1 + c1 * c2 * s3;
            elements(2, 1) = c3 * s2;
            elements(2, 2) = c2;
            elements(2, 3) = s2 * s3;
            elements(3, 1) = -c1 * s3 - c2 * c3 * s1;
            elements(3, 2) = s1 * s2;
            elements(3, 3) = c1 * c3 - c2 * s1 * s3;
        } else if (axis == "ZYZ") {
            elements(1, 1) = c1 * c2 * c3 - s1 * s3;
            elements(1, 2) = -c3 * s1 - c1 * c2 * s3;
            elements(1, 3) = c1 * s2;
            elements(2, 1) = c1 * s3 + c2 * c3 * s1;
            elements(2, 2) = c1 * c3 - c2 * s1 * s3;
            elements(2, 3) = s1 * s2;
            elements(3, 1) = -c3 * s2;
            elements(3, 2) = s2 * s3;
            elements(3, 3) = c2;
        } else if (axis == "ZXZ") {
            elements(1, 1) = c1 * c3 - c2 * s1 * s3;
            elements(1, 1) = -c1 * s3 - c2 * c3 * s1;
            elements(1, 3) = s1 * s2;
            elements(2, 1) = c3 * s1 + c1 * c2 * s3;
            elements(2, 2) = c1 * c2 * c3 - s1 * s3;
            elements(2, 3) = -c1 * s2;
            elements(3, 1) = s2 * s3;
            elements(3, 2) = c3 * s2;
            elements(3, 3) = c2;
        } else {
            throw std::invalid_argument("Unknown axis");
        }

        return linalg::Matrix<T, 3, 3>(elements);
    }

}; // namespace astro 
}; // namespace lao 

#endif // LAO_ASTRO_MATH_TRANSFORMS_H_