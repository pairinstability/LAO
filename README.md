Linear algebra and Astrodynamics Operations library
===

WIP


Linear algebra and other math library for astrodynamics.

Features:
- header-only C++20 linear algebra library.
- supports dense and sparse matrices, which can also be used to construct row and column vectors.
- arithmetic operators are implemented using expression templates.
- matrix size is validated during arithmetic operations at compile-time using concepts.
- various linear algebra solvers.
- planetary positions with JPL low precision.
- date system supporting epochs (JD, MJD, MJD2000).
- various astrodynamics solvers and other utilities.

The build script `scripts/build.sh` can be used to build the examples or unit tests. `scripts/run_tests.sh` runs the unit tests.

Usage, for example:

```cpp
#include <iostream>
#include <lao/lao.hpp>

int main()
{
    lao::linalg::Matrix<double, 1, 2> A { { 5, 6 } };
    lao::linalg::Matrix<double, 1, 2> B { { 2, 3 } };
    lao::linalg::Matrix<double, 1, 2> C { { 3, 1 } };
    lao::linalg::Matrix<double, 1, 2> D;
    lao::linalg::Matrix<double, 1, 2> E;
    D = A + B + C - B;
    E = A - B;

    std::cout << D << std::endl;
    std::cout << E << std::endl;

    lao::linalg::Matrix<double, 2, 3> F { { 1, 2, 1}, {2, 5, 1}};
    lao::linalg::Matrix<double, 3, 2> G { { 5, 6 }, {1, 5}, {2,1}};
    lao::linalg::Matrix<double, 2, 2> H;
    lao::linalg::Matrix<double, 1, 2> I;

    lao::linalg::Matrix<double, 6, 6> temp {{1,2,3,4,5,6}, {7,8,9,10,11,12}, {13, 14, 15, 16, 17, 18}, {19, 20, 21, 22, 23, 24}, {25, 26, 27, 28, 29, 30}, {31, 32, 33, 34, 35, 36}};

    std::cout << *temp.col_begin(4) << std::endl;
    std::cout << *temp.col_end(4) << std::endl;

    std::cout << *temp.row_begin(3) << std::endl;
    std::cout << *temp.row_end(3) << std::endl;

    for (auto x : temp.row_begin(6)) {
        std::cout << "row 6, cols: " << x << std::endl;
    }

    for (auto x : temp.col_begin(3)) {
        std::cout << "col 3, rows: " << x << std::endl;
    }


    H = F * G;
    I = A != B;

    std::cout << H << std::endl;
    std::cout << I << std::endl;

    lao::linalg::Matrix<double, 6, 6>J(lao::linalg::filltype::rand);

   std::cout << J << std::endl;

    auto lambda = []() -> double {
        static double val = 1.0;
        return val++;
    };

    J.fillf(lambda);

    std::cout << J << std::endl;
    std::cout << J(1,2) << std::endl;

    lao::linalg::Matrix<double, 1, 2> K(A);

    std::cout << A << std::endl;
    std::cout << K << std::endl;
    std::cout << 5.0 * K << std::endl;


    lao::linalg::Matrix<double, 3, 3> L {{1,1,2}, {2,1,3}, {3, 1, 1}};
    lao::linalg::Matrix<double, 3, 3> M;
    lao::linalg::Matrix<double, 3, 3> N;

    lao::linalg::LU_doolittle(L, M, N);

    std::cout << "lower: " << std::endl;
    std::cout << M << std::endl;
    std::cout << "upper: " << std::endl;
    std::cout << N << std::endl;

    lao::linalg::Matrix<double, 2, 2> O {{2,1},{5,7}};
    lao::linalg::Matrix<double, 2, 1> P {{11},{13}};
    lao::linalg::Matrix<double, 2, 1> R;

    lao::linalg::solve_jacobi_element(R, O, P, 100, 1e-10);
    std::cout << R << std::endl;

    auto body = lao::astro::JPLLP("Mercury");
    std::cout << "state: " << body.eph(16263.0) << std::endl;
    std::cout << body << std::endl;

    return 0;
}
```