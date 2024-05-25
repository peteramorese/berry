#include "berry/Polynomial.h"
#include "berry/MultiIndex.h"
#include "berry/Operations.h"
#include "berry/BernsteinTransform.h"

#include <Eigen/Dense>

using namespace BRY;

int main() {
    //constexpr std::size_t DIM = 2;
    //Polynomial<DIM> p(2);

    //p.coeff(0, 0) = 1;
    //p.coeff(1, 0) = 2;
    //p.coeff(2, 0) = 3;
    //p.coeff(0, 1) = 4;
    //p.coeff(1, 1) = 5;
    //p.coeff(2, 1) = 6;
    //p.coeff(0, 2) = 7;
    //p.coeff(1, 2) = 8;
    //p.coeff(2, 2) = 9;

    //constexpr std::size_t DIM = 1;
    //Polynomial<DIM> p(2);

    //p.coeff(0) = 1;
    //p.coeff(1) = 1;
    //p.coeff(2) = 1;
    ////p.coeff(3) = 1;

    //auto p_exp = p * p;
    //std::cout << p_exp.tensor() << std::endl;
    //DEBUG(p_exp);
    //DEBUG(p_exp.nMonomials());
    //DEBUG(p_exp.degree());

    constexpr std::size_t DIM = 2;
    Eigen::MatrixXd p_to_b = BernsteinBasisTransform<DIM>::pwrToBernMatrix(7, 1);
    //DEBUG("t mat: \n" << p_to_b);
    Eigen::MatrixXd p_to_b_inc = BernsteinBasisTransform<DIM>::pwrToBernMatrix(3, 5);
    //DEBUG("t mat inc: \n" << p_to_b_inc);

    Eigen::MatrixXd diff =(p_to_b.block(0, 0, p_to_b_inc.rows(), p_to_b_inc.cols()) - p_to_b_inc);
    double d = diff.cwiseAbs().sum();
    DEBUG("Total diff: " << d);

    //bry_int_t test_degree = 2;
    //INFO("Before tmat");
    //Eigen::MatrixXd tmat = BernsteinBasisTransform<DIM>::getTfMatrix(test_degree, 0);
    //INFO("Done!");
    //INFO("T");
    //std::cout << tmat << std::endl;
    //NEW_LINE;

    ////INFO("Before inverse");
    ////Eigen::MatrixXd tmat_inv = tmat.inverse();
    ////INFO("Done!");

    //INFO("Before reverse");
    //Eigen::MatrixXd tmat_rev = BernsteinBasisTransform<DIM>::getInvTfMatrix(test_degree);
    //INFO("Done!");

    //INFO("Before mult");
    //Eigen::MatrixXd tmat_mult = tmat_rev * tmat;
    //INFO("Done!");

    //INFO("I");
    //std::cout << tmat_mult << std::endl;
    //NEW_LINE;

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

    



    //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix(3, 3);
    //matrix << 1, 2, 3,
    //          4, 5, 6,
    //          7, 8, 9;

    //std::cout << matrix << std::endl;

    //// Access the second column (index 1)
    //Eigen::VectorXd secondColumn = matrix.col(1);

    //// Print the second column
    //std::cout << "Second column:\n" << secondColumn << std::endl;

    //Eigen::TensorMap<Eigen::Tensor<double, 1>> t(matrix.col(1).data(), 3);
    //t = 3.0 * t;

    //std::cout<< matrix << std::endl;

}
