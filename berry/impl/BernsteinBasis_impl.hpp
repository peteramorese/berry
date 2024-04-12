#pragma once

#include "BernsteinBasis.h"
#include "MultiIndex.h"
#include "Operations.h"

//template <std::size_t DIM>
//Polynomial<DIM, BRY::Basis::Bernstein> BRY::BernsteinBasis<DIM>::to(const Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_deg_t degree_increase = 0) {
//    
//}

template <std::size_t DIM>
Eigen::MatrixXd BRY::BernsteinBasis<DIM>::getTransformationMatrix(const Polynomial<DIM, Basis::Power>& p, bry_deg_t degree_increase) {
    bry_deg_t to_degree = p.degree() + degree_increase;
    Eigen::MatrixXd matrix(p.nMonomials(), p.nMonomials());
    matrix.setZero();
    
    // In place multi index for 'l' indices
    std::array<bry_idx_t, DIM> l_idx;
    std::array<bry_idx_t, DIM> i_idx;

    //// Compute the denominator coefficients
    //std::vector<bry_idx_t> denom_coeffs;
    //denom_coeffs.reserve(p.nMonomials());
    //for (auto midx = mIdx(DIM, p.degree() + 1); !midx.last(); ++midx) {
    //    denom_coeffs.push_back();
    //}

    MultiIndex<ExhaustiveIncrementerWrap> i_midx(i_idx.data(), DIM, true, p.degree() + 1);
    for (; !i_midx.last(); ++i_midx) {
        
        // Use the row midx as the bounds for the column (i) iterator
        std::vector<bry_idx_t> index_bounds;
        index_bounds.reserve(DIM);
        for (bry_idx_t r_i : i_midx)
            index_bounds.push_back(r_i + 1);

        DEBUG("imidx: " << i_midx);

        // Create the column multi index
        MultiIndex<BoundedExhaustiveIncrementerWrap> l_midx(l_idx.data(), DIM, true, index_bounds, p.degree() + 1);
        for (; !l_midx.last(); ++l_midx) {
            
            DEBUG("     lmidx: " << l_midx);

            // Compute the transformation coefficient in a numerically stable way
            bry_float_t transformation_coeff = 1.0;
            for (bry_idx_t j = 0; j < DIM; ++j) {
                DEBUG("         denominator input: " << to_degree << " bottom : " << l_midx[j]);
                DEBUG("         denominator: " << static_cast<bry_float_t>(binom(to_degree, l_midx[j])));
                transformation_coeff *= static_cast<bry_float_t>(binom(i_midx[j], l_midx[j])) / static_cast<bry_float_t>(binom(to_degree, l_midx[j]));
            }

            matrix(i_midx.inc().wrappedIdx(), l_midx.inc().wrappedIdx()) = transformation_coeff;
        }
    }
    return matrix;
}

//Eigen::Tensor<bry_float_t, DIM> getDenominatorTensor()