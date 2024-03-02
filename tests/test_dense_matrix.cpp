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

/// @brief Test default constructor.
TEST_F(MatrixTest, DefaultConstructor)
{
    lao::linalg::Matrix<int, 50, 50> mat;
    lao::linalg::Matrix<int, 50, 50> mat_zeros;
    mat_zeros.fill(0);

    EXPECT_EQ(true, matricesEqual(mat, mat_zeros));
}

/// @brief Test passing a matrix to the constructor.
TEST_F(MatrixTest, MatrixConstructor)
{
    lao::linalg::Matrix<int, 2, 2> mat({ { 1, 2 }, { 3, 4 } });
    EXPECT_EQ(mat(1, 1), 1);
    EXPECT_EQ(mat(1, 2), 2);
    EXPECT_EQ(mat(2, 1), 3);
    EXPECT_EQ(mat(2, 2), 4);
}

/// @brief Test fill with zeros in constructor.
TEST_F(MatrixTest, FillZeros)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::zeros);

    for (size_t i = 1; i <= mat.rows(); ++i) {
        for (size_t j = 1; j <= mat.cols(); ++j) {
            EXPECT_EQ(mat(i, j), 0);
        }
    }
}

/// @brief Test fill with ones in constructor.
TEST_F(MatrixTest, FillOnes)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::ones);

    for (size_t i = 1; i <= mat.rows(); ++i) {
        for (size_t j = 1; j <= mat.cols(); ++j) {
            EXPECT_EQ(mat(i, j), 1);
        }
    }
}

/// @brief Test fill with eye.
TEST_F(MatrixTest, FillEye)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::eye);

    for (size_t i = 1; i <= mat.rows(); ++i) {
        for (size_t j = 1; j <= mat.cols(); ++j) {
            if (i == j) {
                EXPECT_EQ(mat(i, j), 1);
            } else {
                EXPECT_EQ(mat(i, j), 0);
            }
        }
    }
}

/// @brief Test copy constructor.
TEST_F(MatrixTest, CopyConstructor)
{
    lao::linalg::Matrix<int, 3, 3> mat1;
    mat1.fill(5);
    lao::linalg::Matrix<int, 3, 3> mat2(mat1);

    EXPECT_TRUE(matricesEqual(mat1, mat2));
}

/// @brief Test copy assignment operator.
TEST_F(MatrixTest, CopyAssignmentOperator)
{
    lao::linalg::Matrix<int, 3, 3> mat1;
    mat1.fill(5);
    lao::linalg::Matrix<int, 3, 3> mat2;
    mat2 = mat1;

    EXPECT_TRUE(matricesEqual(mat1, mat2));
}

/// @brief Test row iterators.
TEST_F(MatrixTest, RowIterators)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::ones);

    for (size_t i = 1; i <= mat.rows(); ++i) {
        int row_sum = 0;
        for (auto it = mat.row_begin(i); it != mat.row_end(i); ++it) {
            row_sum += *it;
        }
        EXPECT_EQ(row_sum + 1, mat.cols());
    }
}

/// @brief Test column iterators.
TEST_F(MatrixTest, ColumnIterators)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::ones);

    for (size_t j = 1; j <= mat.cols(); ++j) {
        int col_sum = 0;
        for (auto it = mat.col_begin(j); it != mat.col_end(j); ++it) {
            col_sum += *it;
        }
        EXPECT_EQ(col_sum + 1, mat.rows());
    }
}

/// @brief Test fillf.
TEST_F(MatrixTest, FillfMethod)
{
    lao::linalg::Matrix<int, 3, 3> mat;
    mat.fillf([]() { return 5; });

    for (size_t i = 1; i <= mat.rows(); ++i) {
        for (size_t j = 1; j <= mat.cols(); ++j) {
            EXPECT_EQ(mat(i, j), 5);
        }
    }
}

/// @brief Test reset.
TEST_F(MatrixTest, ResetMethod)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::ones);

    EXPECT_FALSE(mat.is_empty());

    mat.reset();

    EXPECT_TRUE(mat.is_empty());
}

/// @brief Test accessing an element and modifying.
TEST_F(MatrixTest, ElementAccessAndModification)
{
    lao::linalg::Matrix<int, 3, 3> mat(lao::linalg::filltype::zeros);

    for (size_t i = 1; i <= mat.rows(); ++i) {
        for (size_t j = 1; j <= mat.cols(); ++j) {
            mat(i, j) = i + j;
            EXPECT_EQ(mat(i, j), i + j);
        }
    }
}

/// @brief Test accessing elements from iterators.
TEST_F(MatrixTest, ElementAccessFromIterator)
{
    lao::linalg::Matrix<double, 6, 6> temp { { 1, 2, 3, 4, 5, 6 }, { 7, 8, 9, 10, 11, 12 }, { 13, 14, 15, 16, 17, 18 }, { 19, 20, 21, 22, 23, 24 }, { 25, 26, 27, 28, 29, 30 }, { 31, 32, 33, 34, 35, 36 } };

    EXPECT_EQ(*temp.col_begin(4), 4);
    EXPECT_EQ(*temp.col_end(4), 34);

    EXPECT_EQ(*temp.row_begin(3), 13);
    EXPECT_EQ(*temp.row_end(3), 18);

    std::vector<int> scratch_row6;
    std::vector<int> scratch_col2;
    std::vector<int> expect_row6 = { 31, 32, 33, 34, 35, 36 };
    std::vector<int> expect_col2 = { 2, 8, 14, 20, 26, 32 };

    for (auto x : temp.row_begin(6)) {
        scratch_row6.push_back(x);
    }

    for (auto x : temp.col_begin(2)) {
        scratch_col2.push_back(x);
    }

    EXPECT_TRUE(std::equal(scratch_row6.begin(), scratch_row6.end(), expect_row6.begin()));
    EXPECT_TRUE(std::equal(scratch_col2.begin(), scratch_col2.end(), expect_col2.begin()));
}

/// @brief Test accessing elements out of bounds.
TEST_F(MatrixTest, AccessOutOfBounds)
{
    lao::linalg::Matrix<int, 3, 3> mat;

    EXPECT_THROW(mat(0, 0), std::out_of_range);
    EXPECT_THROW(mat(4, 2), std::out_of_range);
    EXPECT_THROW(mat(1, 4), std::out_of_range);
}

/// @brief Test accessing row iterators out of bounds.
TEST_F(MatrixTest, RowIteratorsOutOfRange)
{
    lao::linalg::Matrix<int, 3, 3> mat;

    EXPECT_THROW(mat.row_begin(0), std::out_of_range);
    EXPECT_THROW(mat.row_begin(4), std::out_of_range);

    EXPECT_THROW(mat.row_end(0), std::out_of_range);
    EXPECT_THROW(mat.row_end(4), std::out_of_range);
}

/// @brief Test accessing column iterators out of bounds.
TEST_F(MatrixTest, ColIteratorsOutOfRange)
{
    lao::linalg::Matrix<int, 3, 3> mat;

    EXPECT_THROW(mat.col_begin(0), std::out_of_range);
    EXPECT_THROW(mat.col_begin(4), std::out_of_range);

    EXPECT_THROW(mat.col_end(0), std::out_of_range);
    EXPECT_THROW(mat.col_end(4), std::out_of_range);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
