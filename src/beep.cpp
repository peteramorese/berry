#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"
#include "BernsteinBasis.h"

using namespace BRY;

int main() {
    //Polynomial<3> p1(3);
    Polynomial<1> p(3);

    p.coeff(0) = 1;
    p.coeff(1) = 2;
    p.coeff(2) = 3;
    p.coeff(3) = 4;

    DEBUG(p);

    std::cout << BernsteinBasis<1>::getTransformationMatrix(p) << std::endl;

    for (int i = 0; i < 5; ++ i) {
        auto pr = pascalRow(i);
        //INFO("Pascals row " << i);
        for (auto n : pr) 
            std::cout << n << " ";
        std::cout << std::endl;
    }
}
