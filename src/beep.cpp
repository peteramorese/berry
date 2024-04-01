#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"

using namespace BRY;

int main() {
    Polynomial<3> p1(3);
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

    //p1.coeff(1, 3, 2) = 2.0;
    //p1.coeff(0, 0, 0) = 3.0;
    //p1.coeff(2, 0, 3) = 4.0;

    p1.coeff(1, 0, 0) = 1.0;
    p1.coeff(0, 1, 0) = 1.0;
    p1.coeff(0, 0, 1) = 1.0;

    p2.coeff(1, 3, 2) = 5.0;
    p2.coeff(0, 2, 3) = 4.0;
    p2.coeff(2, 4, 5) = 3.0;

    DEBUG("p1: " << p1);
    DEBUG("p2: " << p2);
    auto p_sum = 66.0 + p1 + p2;
    auto p_diff = p2 - 2.5 * p1;
    auto p_mult = p1 * p2;
    auto p_exp = p1 ^ 3;
    //DEBUG("p1 + p2: " << p_sum);
    //DEBUG("p1 - p2: " << p_diff);
    DEBUG("p1 * p2: " << p_mult);
    DEBUG("p1 ^ 3: " << p_exp);
    //DEBUG("p2 ^ 3: " << p2 ^ 5);

    //Eigen::MatrixXd m(2, 2);
    //m(0, 0) = 1;
    //m(1, 0) = 1;
    //m(0, 1) = 1;
    //m(1, 1) = 1;

    //std::cout << m << std::endl;
    //m.conservativeResize(3,3);

    //std::cout << m << std::endl;

    Eigen::Tensor<bry_float_t, 2> t(3, 3);
    t.setValues({{1,2,3}, {4,5,6}, {7,8,9}});
    std::cout << t <<std::endl;
    std::array<int64_t, 2> offsets = {{1, 1}};
    std::array<int64_t, 2> extents = {{2, 2}};
    //Eigen::Tensor<bry_float_t, 2>
    std::cout << t.slice(offsets, extents) <<std::endl;
}
