#include "Polynomial.h"

using namespace BRY;

int main() {
    Polynomial<3> p1(4);
    Polynomial<4> p2(4);

    std::array<bry_deg_t, 3> p1_exps;
    p1_exps[0] = 1;
    p1_exps[1] = 2;
    p1_exps[2] = 0;

    std::array<bry_deg_t, 4> p2_exps;
    p2_exps[0] = 1;
    p2_exps[1] = 2;
    p2_exps[2] = 0;
    p2_exps[3] = 0;

    DEBUG(p1.wrap(p1_exps));
    DEBUG(p2.wrap(p2_exps));
    p1.coeff(1, 3, 3) = 54.0;
    p1.coeff(0, 0, 1) = 4.0;
    DEBUG(p1.coeff(1, 3, 3));

    DEBUG(p1);
    PRINT_NAMED("beep", "borp");
}