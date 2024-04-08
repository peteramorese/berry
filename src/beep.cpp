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

    // Create 2 matrices using tensors of rank 2
    //Eigen::Tensor<int, 3> a(2, 2, 2);
    //a.setValues({{{1, 5}, {3, 7}}, {{2, 6}, {4, 8}}});
    //Eigen::Tensor<int, 2> a_vec = a.reshape(std::array<int, 2>{{8, 1}});
    ////Eigen::Tensor<int, 1> b_mat = b.reshape(std::array<int, 2>{{4, 4}})

    //for (MultiIndex midx(3, 2); !midx.right(); ++midx) {
    //    std::array<int, 3> inp;
    //    for (std::size_t i = 0; i < 3; ++i) {
    //        inp[i] = midx[i];
    //    }
    //    DEBUG("Midx: " << midx << " element: " << a(inp));
    //}

    //std::cout << a << std::endl;
    //std::cout << a_vec << std::endl;
    //std::cout << b << std::endl;

    //// Compute the traditional matrix product
    //Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>(1, 0) };
    //Eigen::Tensor<int, 2> prod = a.contract(b, product_dims);

    //std::cout << prod << std::endl;
    //BoundedExhaustiveIncrementer incrementer(std::vector<bry_idx_t>{2, 1, 4}, 4);
    for (auto midx = mIdxBEW(std::vector<bry_idx_t>{2, 1, 4}, 4); !midx.last(); ++midx) {
        DEBUG("Midx: " << midx << " wrapped idx: " << midx.incrementer().wrappedIdx());
    }

    uint32_t i = 0;
    for (auto midx = mIdxW(3, 4); !midx.last(); ++midx) {
        DEBUG("Exhaustive Midx: " << midx << " wrapped idx: " << i++ << " comp to: " << midx.incrementer().wrappedIdx());
    }

}
