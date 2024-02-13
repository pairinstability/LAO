LAO - Linear Algebra Operations
===

Note this isn't industrial-grade, this is just a fun project.

Features:
- header-only C++20 linear algebra library.
- supports dense and sparse matrices, which can also be used to construct row and column vectors.
- arithmetic operators are implemented using expression templates.
- matrix size is validated during arithmetic operations at compile-time using concepts.

The build script `scripts/build.sh` can be used to build the examples or unit tests. `scripts/run_tests.sh` runs the unit tests.

Usage, for example:

```cpp
#include <iostream>
#include <lao/lao.hpp>

int main()
{
    lao::Matrix<double, 1, 2> A { { 5, 6 } };
    lao::Matrix<double, 1, 2> B { { 2, 3 } };
    lao::Matrix<double, 1, 2> C { { 3, 1 } };
    lao::Matrix<double, 1, 2> D;
    lao::Matrix<double, 1, 2> E;
    D = A + B + C - B;
    E = A - B;

    std::cout << D << std::endl;
    std::cout << E << std::endl;

    lao::Matrix<double, 2, 3> F { { 1, 2, 1}, {2, 2, 1}};
    lao::Matrix<double, 3, 2> G { { 5, 6 }, {1, 5}, {2,1}};
    lao::Matrix<double, 2, 2> H;
    lao::Matrix<double, 1, 2> I;

    H = F * G;
    I = A != B;

    std::cout << H << std::endl;
    std::cout << I << std::endl;

    lao::Matrix<double, 6, 6>J(lao::filltype::rand);

    std::cout << J << std::endl;

    auto lambda = []() -> double {
        static double val = 1.0;
        return val++;
    };

    J.fillf(lambda);

    std::cout << J << std::endl;
    std::cout << J(1,2) << std::endl;

    lao::Matrix<double, 1, 2> K(A);

    std::cout << A << std::endl;
    std::cout << K << std::endl;

    return 0;
}
```