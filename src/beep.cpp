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

    //BRY::Polynomial<3> p(2);
    //p.coeff(0, 0, 0) = 100;
    //p.coeff(1, 0, 0) = 1;
    //p.coeff(0, 1, 0) = 55;
    //p.coeff(1, 1, 0) = 2;
    //p.coeff(1, 2, 0) = 3;
    //p.coeff(0, 2, 1) = 4;
    //p.coeff(0, 1, 1) = 5;
    //p.coeff(0, 1, 2) = 6;
    //p.coeff(0, 0, 1) = 7;
    //p.coeff(2, 0, 1) = 8;
    //p.coeff(2, 2, 2) = 9;
    //bry_float_t result = p(.3, .4, .5);
    //INFO(p(.3, .4, .5));
    //INFO(p(0, 0, 0));
    ////std::cout << std::fixed << std::setprecision(10) << result << std::endl;

    //INFO("p: " << p);
    //Polynomial<3> p_prime = p.derivative(2);
    //INFO("p derivative: " << p_prime);

    //BRY::Polynomial<2> p1(0);
    //p1.coeff(0, 0) = 1;
    //BRY::Polynomial<2> p2(2);
    //p2.coeff(0, 0) = 10;
    //p2.coeff(1, 0) = 1;
    //p2.coeff(0, 1) = 6;
    //DEBUG("p1: " << p1);
    //DEBUG("p2: " << p2);
    //auto p3 = p1 * p2;
    //DEBUG("p3: " << p3);

    //BRY::Polynomial<1> f(2);
    //f.coeff(0) = 2;
    //f.coeff(1) = 3;
    //f.coeff(2) = 4;
    //INFO("f^0: " << (f^0));
    //INFO("f^1: " << (f^1));
    //INFO("f^2: " << (f^2));
    //NEW_LINE;

    BRY::Polynomial<2> p(2);
    p.coeff(0, 0) = 2;
    p.coeff(1, 0) = 3;
    p.coeff(2, 0) = 4;
    p.coeff(0, 1) = 1;
    DEBUG("p: " << p);
    auto p_lifted = p.liftDegree(4);
    DEBUG("p raised deg: " << p_lifted << " deg: " << p_lifted.degree());
    //INFO("(f + v)^0: " << (p^0));
    //INFO("(f + v)^1: " << (p^1));
    //INFO("(f + v)^2: " << (p^2));

    Matrix tf = makeDegreeChangeTransform<2>(2, 4);
    DEBUG("tf: \n" << tf);
    auto p_lifted_2 = transform(p, tf);
    DEBUG("p raised deg 2: " << p_lifted_2 << " deg: " << p_lifted_2.degree());
}
