#include <lao/lao.hpp>
#include <iostream>

int main()
{
    lao::Matrix<double, 2, 2> a{ {1, 2}, {2,2} };
    lao::Matrix<double, 3, 3> b{ {9, 8, 3}, {2,1,4}, {4,0,9}};
    lao::Matrix<double, 3, 3> c(lao::filltype::rand);

    std::cout << a(1,1) << std::endl;
    std::cout << b << std::endl;
    std::cout << c << std::endl;


    auto lambda = []() -> double {
        static double val = 1.0;
        return val++;
    };

    c.fillf(lambda);

    std::cout << c << std::endl;



    return 0;
}