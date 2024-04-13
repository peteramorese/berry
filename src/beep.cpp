#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"
#include "BernsteinBasis.h"

#include <Eigen/Dense>

using namespace BRY;

int main() {
    constexpr std::size_t DIM = 2;
    Polynomial<DIM> p(2);

    p.coeff(0, 0) = 1;
    p.coeff(1, 0) = 2;
    p.coeff(2, 0) = 3;
    p.coeff(0, 1) = 4;
    p.coeff(1, 1) = 5;
    p.coeff(2, 1) = 6;
    p.coeff(0, 2) = 7;
    p.coeff(1, 2) = 8;
    p.coeff(2, 2) = 9;

    //constexpr std::size_t DIM = 1;
    //Polynomial<DIM> p(3);

    //p.coeff(0) = 1;
    //p.coeff(1) = 2;
    //p.coeff(2) = 3;
    //p.coeff(3) = 4;

    DEBUG(p);
    DEBUG(p.nMonomials());
    DEBUG(p.degree());


    Eigen::MatrixXd tmat = BernsteinBasis<DIM>::getTransformationMatrix(p.degree());
    NEW_LINE;
    INFO("T");
    std::cout << tmat << std::endl;
    NEW_LINE;
    INFO("T^-1");
    //std::cout << tmat.inverse().cwiseAbs() << std::endl;
    std::cout << tmat.inverse() << std::endl;
    NEW_LINE;
    INFO("Backwards");
    //std::cout << BernsteinBasis<DIM>::getInverseTransformationMatrix(p.degree()).cwiseAbs() << std::endl;
    std::cout << BernsteinBasis<DIM>::getInverseTransformationMatrix(p.degree()) << std::endl;

}
