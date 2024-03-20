#include "Polynomial.h"

using namespace BRY;

int main() {
    Polynomial<3> poly(4);

    poly.coeff(1, 3, 3) = 54.0;
    DEBUG(poly.coeff(1,3,3));

}