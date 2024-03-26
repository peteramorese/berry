#include "Polynomial.h"

using namespace BRY;

int main() {
    Polynomial<3> p1(5);
    Polynomial<3> p2(5);

    //std::array<bry_deg_t, 2> p1_exps;
    //p1_exps[0] = 1;
    //p1_exps[1] = 2;
    ////p1_exps[2] = 0;

    //std::array<bry_deg_t, 2> p2_exps;
    //p2_exps[0] = 1;
    //p2_exps[1] = 2;
    ////p2_exps[2] = 0;

    //DEBUG(p1.wrap(p1_exps));
    //DEBUG(p2.wrap(p1_exps));

    p1.coeff(1, 3, 2) = 54.0;
    p1.coeff(0, 0, 0) = 4.0;
    p1.coeff(2, 0, 3) = 9.0;

    p2.coeff(1, 3, 2) = 46.0;
    p2.coeff(0, 2, 3) = 9.0;
    p2.coeff(2, 4, 5) = 9.0;

    DEBUG("p1: " << p1);
    DEBUG("p2: " << p2);
    auto p_sum = p1 + 5.0 * p2;
    auto p_diff = p1 - p2;
    auto p_mult = p1 * p2;
    DEBUG("p1 + p2: " << p_sum);
    DEBUG("p1 - p2: " << p_diff);
    DEBUG("p1 * p2: " << p_mult);
}
