#include <cstddef>
#include <gtest/gtest.h>
#include <lao/lao.hpp>
#include <sstream>
#include <string>

class AstroTest : public ::testing::Test {
protected:
    template <typename S, size_t R, size_t C>
    bool matricesEqual(const lao::linalg::Matrix<S, R, C>& mat1, const lao::linalg::Matrix<S, R, C>& mat2) const
    {
        double tol = 1e-9;
        for (size_t i = 1; i <= mat1.rows(); ++i) {
            for (size_t j = 1; j <= mat1.cols(); ++j) {
                if (fabs(mat1(i,j) - mat2(i,j)) < tol)
                    return false;
            }
        }
        return true;
    }
};

TEST_F(AstroTest, JPLLP_Ephemeris)
{
    lao::linalg::Matrix<double, 1, 6> expected_state {{3.92554e+10, -1.47559e+10, -2.27302e+10, 21964.7, 13898.2, 39307.7}};
    auto body = lao::astro::JPLLP("Mercury");

    EXPECT_EQ(matricesEqual(body.eph(16263.0), expected_state), true);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
