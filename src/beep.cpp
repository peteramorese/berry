#include "berry/Polynomial.h"
#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "berry/BernsteinTransform.h"

#include <Eigen/Dense>

#include <iomanip>

using namespace BRY;

int main() {

    //BRY::Polynomial<2> p(1);
    //p.coeff(0, 0) = 100;
    //p.coeff(1, 0) = 2;
    //p.coeff(0, 1) = 3;
    //p.coeff(1, 1) = 4;
    //bry_float_t result = p(.3, .4);
    //INFO("I got: " << result);
    //INFO("Should be: 102.28");

    BRY::Polynomial<3> p(2);
    p.coeff(0, 0, 0) = 100;
    p.coeff(1, 0, 0) = 1;
    p.coeff(0, 1, 0) = 55;
    p.coeff(1, 1, 0) = 2;
    p.coeff(1, 2, 0) = 3;
    p.coeff(0, 2, 1) = 4;
    p.coeff(0, 1, 1) = 5;
    p.coeff(0, 1, 2) = 6;
    p.coeff(0, 0, 1) = 7;
    p.coeff(2, 0, 1) = 8;
    p.coeff(2, 2, 2) = 9;
    bry_float_t result = p(.3, .4, .5);
    INFO(p(.3, .4, .5));
    INFO(p(0, 0, 0));
    //std::cout << std::fixed << std::setprecision(10) << result << std::endl;

    INFO("p: " << p);
    Polynomial<3> p_prime = p.derivative(2);
    INFO("p derivative: " << p_prime);
}
