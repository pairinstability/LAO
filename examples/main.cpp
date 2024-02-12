#include <iostream>
#include <lao/lao.hpp>

int main()
{
    lao::Matrix<double, 1, 2> A { { 1, 2 } };
    lao::Matrix<double, 1, 2> B { { 2, 3 } };
    lao::Matrix<double, 1, 2> C { { 3, 4 } };
    lao::Matrix<double, 1, 2> S;
    S = A + B + C;

    std::cout << S << std::endl;

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