#include <iostream>
#include <lao/lao.hpp>

int main()
{
    lao::Matrix<double, 1, 2> A { { 5, 6 } };
    lao::Matrix<double, 1, 2> B { { 2, 3 } };
    lao::Matrix<double, 1, 2> C { { 3, 1 } };
    lao::Matrix<double, 1, 2> S;
    lao::Matrix<double, 1, 2> E;
    S = A + B + C - B;
    E = A - B;

    std::cout << S << std::endl;
    std::cout << E << std::endl;

    /*
        auto lambda = []() -> double {
            static double val = 1.0;
            return val++;
        };

        c.fillf(lambda);

    //    std::cout << c << std::endl;
    */

    return 0;
}