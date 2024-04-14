#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"
#include "BernsteinTransform.h"

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


    bry_deg_t test_degree = 2;
    INFO("Before tmat");
    Eigen::MatrixXd tmat = BernsteinBasisTransform<DIM>::getTfMatrix(test_degree, 0);
    INFO("Done!");
    INFO("T");
    std::cout << tmat << std::endl;
    NEW_LINE;

    //INFO("Before inverse");
    //Eigen::MatrixXd tmat_inv = tmat.inverse();
    //INFO("Done!");

    INFO("Before reverse");
    Eigen::MatrixXd tmat_rev = BernsteinBasisTransform<DIM>::getInvTfMatrix(test_degree);
    INFO("Done!");

    INFO("Before mult");
    Eigen::MatrixXd tmat_mult = tmat_rev * tmat;
    INFO("Done!");

    INFO("I");
    std::cout << tmat_mult << std::endl;
    NEW_LINE;

    //NEW_LINE;
    //INFO("T");
    //std::cout << tmat << std::endl;
    //NEW_LINE;
    //INFO("T^-1");
    ////std::cout << tmat.inverse().cwiseAbs() << std::endl;
    //std::cout << tmat.inverse() << std::endl;
    //NEW_LINE;
    //INFO("Backwards");
    ////std::cout << BernsteinBasis<DIM>::getInverseTransformationMatrix(p.degree()).cwiseAbs() << std::endl;
    //std::cout << BernsteinBasis<DIM>::getInverseTransformationMatrix(p.degree()) << std::endl;

}
