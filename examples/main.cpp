#include "lao/dense/LU.hpp"
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
    std::cout << 5.0 * K << std::endl;


    lao::Matrix<double, 3, 3> L {{1,1,2}, {2,1,3}, {3, 1, 1}};
    lao::Matrix<double, 3, 3> M;
    lao::Matrix<double, 3, 3> N;

    lao::LU_doolittle(L, M, N);

    std::cout << M << std::endl;
    std::cout << N << std::endl;


    return 0;
}