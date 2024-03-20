#include "Polynomial.h"

using namespace BRY;

int main() {
    Polynomial<3> poly(4);

    std::array<bry_deg_t, 3> degrees;
    degrees[0] = 3;
    degrees[1] = 0;
    degrees[2] = 1;
    DEBUG(poly.wrap(degrees));

    std::array<bry_deg_t, 3> got_degrees = poly.unwrap(poly.wrap(degrees));
    PRINT_VEC3("hello ", got_degrees);
}