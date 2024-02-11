/// test_dense_matrix.cpp tests the dense matrix class, and its fundamental operations.
/// This excludes vector specializations and all mathematical operations.

#include <gtest/gtest.h>
#include <lao/lao.hpp>
#include <cstddef>
#include <string>

template <typename MatrixType>
class DenseMatrixTest : public ::testing::Test {
protected:
    void SetUp() override {
        matrix_ones.ones();
    }

    void TearDown() override {
        matrix_ones.reset();
    }

    std::string matrixToString(const MatrixType& matrix) {
        std::ostringstream oss;
        matrix.print(oss);
        return oss.str();
    }

    bool matricesEqual(const MatrixType& mat1, const MatrixType& mat2) const {
        if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols()) {
            return false;
        }
        for (size_t i = 0; i < mat1.rows(); ++i) {
            for (size_t j = 0; j < mat1.cols(); ++j) {
                if (mat1(i, j) != mat2(i, j)) {
                    return false;
                }
            }
        }
        return true;
    }

    MatrixType matrix_zeros;
    MatrixType matrix_ones;
    MatrixType matrix_rand;
};

using DenseMatrixTypes = ::testing::Types<lao::Matrix<double, 50, 50>, lao::Matrix<double, 10, 2>, lao::Matrix<double, 2, 3>>;
TYPED_TEST_SUITE(DenseMatrixTest, DenseMatrixTypes);

/// @brief Test default constructor.
/// @details The default constructor zero-initializes a matrix.
TYPED_TEST(DenseMatrixTest, DefaultConstructor) {
    TypeParam matrix;
    for (size_t i = 0; i < matrix.rows(); ++i) {
        for (size_t j = 0; j < matrix.cols(); ++j) {
            EXPECT_EQ(matrix(i, j), 0.0);
        }
    }
}

/// @brief Test passing a vector to the constructor.
TEST(DenseMatrixTest, VectorConstructor) {
    lao::Matrix<double, 2, 2> mat({{1,2},{3,4}});
    EXPECT_EQ(mat.shape(), 4);
    EXPECT_EQ(mat(0,0), 1);
    EXPECT_EQ(mat(0,1), 2);
    EXPECT_EQ(mat(1,0), 3);
    EXPECT_EQ(mat(1,1), 4);
}


/// @brief Test initializing a matrix with random values.
TYPED_TEST(DenseMatrixTest, InitializationWithRandomValues) {
    TypeParam matrix(lao::filltype::rand);
    EXPECT_FALSE(matrix.is_empty());
}


/// @brief Test using zeros method.
TYPED_TEST(DenseMatrixTest, ZerosMethod) {
    TypeParam matrix;
    matrix.zeros();
    for (size_t i = 0; i < matrix.rows(); ++i) {
        for (size_t j = 0; j < matrix.cols(); ++j) {
            EXPECT_EQ(matrix(i, j), 0.0);
        }
    }
}

/// @brief Test ones method.
TYPED_TEST(DenseMatrixTest, OnesMethod) {
    TypeParam matrix;
    matrix.ones();
    for (size_t i = 0; i < matrix.rows(); ++i) {
        for (size_t j = 0; j < matrix.cols(); ++j) {
            EXPECT_EQ(matrix(i, j), 1);
        }
    }
}

/// @brief Test copying a matrix.
TYPED_TEST(DenseMatrixTest, CopyConstructorAndAssignment) {
    TypeParam matrix_original(lao::filltype::rand);
    // copy constructor.
    TypeParam matrix_copy = matrix_original;
    TypeParam matrix_assigned;
    // copy assignment operator.
    matrix_assigned = matrix_original;
    EXPECT_TRUE(this->matricesEqual(matrix_original, matrix_copy));
    EXPECT_TRUE(this->matricesEqual(matrix_original, matrix_assigned));
}

/// @brief Test accessing an element and modifying.
TYPED_TEST(DenseMatrixTest, ElementAccessAndModification) {
    TypeParam matrix;
    matrix.zeros();
    for (size_t i = 0; i < matrix.rows(); ++i) {
        for (size_t j = 0; j < matrix.cols(); ++j) {
            matrix(i, j) = i + j;
            EXPECT_EQ(matrix(i, j), i + j);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
