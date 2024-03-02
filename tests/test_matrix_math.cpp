#include <cstddef>
#include <gtest/gtest.h>
#include <lao/lao.hpp>
#include <string>

class MatrixTest : public ::testing::Test {
protected:
    template <typename S, size_t R, size_t C>
    bool matricesEqual(const lao::linalg::Matrix<S, R, C>& mat1, const lao::linalg::Matrix<S, R, C>& mat2) const
    {
        for (size_t i = 1; i <= mat1.rows(); ++i) {
            for (size_t j = 1; j <= mat1.cols(); ++j) {
                if (mat1(i, j) != mat2(i, j))
                    return false;
            }
        }
        return true;
    }
};

/// @brief Test matrix addition.
TEST_F(MatrixTest, MatrixAddition)
{
    lao::linalg::Matrix<int, 2, 2> mat1 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat2 { { 5, 6 }, { 7, 8 } };
    lao::linalg::Matrix<int, 2, 2> mat3;

    lao::linalg::Matrix<int, 2, 2> result { { 6, 8 }, { 10, 12 } };

    mat3 = mat1 + mat2;

    EXPECT_EQ(matricesEqual(mat3, result), true);
}

/// @brief Test matrix subtraction.
TEST_F(MatrixTest, MatrixSubtraction)
{
    lao::linalg::Matrix<int, 2, 2> mat1 { { 5, 6 }, { 7, 8 } };
    lao::linalg::Matrix<int, 2, 2> mat2 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat3;

    lao::linalg::Matrix<int, 2, 2> result { { 4, 4 }, { 4, 4 } };

    mat3 = mat1 - mat2;

    EXPECT_EQ(matricesEqual(mat3, result), true);
}

/// @brief Test matrix multiplication.
TEST_F(MatrixTest, MatrixMultiplication)
{
    lao::linalg::Matrix<int, 2, 2> mat1 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat2 { { 5, 6 }, { 7, 8 } };
    lao::linalg::Matrix<int, 2, 2> mat3;

    lao::linalg::Matrix<int, 2, 2> result { { 19, 22 }, { 43, 50 } };

    mat3 = mat1 * mat2;

    EXPECT_EQ(matricesEqual(mat3, result), true);
}

/// @brief Test combination of matrix operations.
TEST_F(MatrixTest, CombinedOperations)
{
    lao::linalg::Matrix<int, 2, 2> mat1 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat2 { { 5, 6 }, { 7, 8 } };
    lao::linalg::Matrix<int, 2, 2> mat3 { { 9, 10 }, { 11, 12 } };
    lao::linalg::Matrix<int, 2, 2> mat4;

    lao::linalg::Matrix<int, 2, 2> result { { 418, 460 }, { 944, 1038 } };

    mat4 = mat1 * mat2 * mat3 + mat2;

    EXPECT_EQ(matricesEqual(mat4, result), true);
}

/// @brief Test scalar multiplication.
TEST_F(MatrixTest, ScalarMultiplication)
{
    lao::linalg::Matrix<int, 2, 2> mat { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat2;

    lao::linalg::Matrix<int, 2, 2> result { { 2, 4 }, { 6, 8 } };

    mat2 = 2 * mat;

    EXPECT_EQ(matricesEqual(mat2, result), true);
}

/// @brief Test element-wise matrix multiplication.
TEST_F(MatrixTest, ElementWiseMultiplication)
{
    lao::linalg::Matrix<int, 2, 2> mat1 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat2 { { 5, 6 }, { 7, 8 } };
    lao::linalg::Matrix<int, 2, 2> mat3;

    lao::linalg::Matrix<int, 2, 2> result { { 5, 12 }, { 21, 32 } };

    mat3 = mat1 % mat2;

    EXPECT_EQ(matricesEqual(mat3, result), true);
}

/// @brief Test element-wise equality.
TEST_F(MatrixTest, ElementWiseEquality)
{
    lao::linalg::Matrix<int, 2, 2> mat1 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat2 { { 1, 2 }, { 3, 4 } };
    lao::linalg::Matrix<int, 2, 2> mat3 { { 5, 6 }, { 7, 8 } };

    EXPECT_EQ(matricesEqual(mat1, mat2), true);
    EXPECT_EQ(matricesEqual(mat1, mat3), false);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
