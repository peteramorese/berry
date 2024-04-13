#pragma once

#include "BernsteinBasis.h"
#include "MultiIndex.h"
#include "Operations.h"

//template <std::size_t DIM>
//Polynomial<DIM, BRY::Basis::Bernstein> BRY::BernsteinBasis<DIM>::to(const Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_deg_t degree_increase = 0) {
//    
//}

template <std::size_t DIM>
Eigen::MatrixXd BRY::BernsteinBasis<DIM>::getTransformationMatrix(bry_deg_t degree, bry_deg_t degree_increase) {
    bry_deg_t to_degree = degree + degree_increase;

    auto makeCoeff = [&] (const auto& i_midx, const auto& l_midx) -> bry_float_t {
        // Compute the transformation coefficient in a numerically stable way
        bry_float_t transformation_coeff = 1.0;
        for (bry_idx_t j = 0; j < DIM; ++j) {
            transformation_coeff *= static_cast<bry_float_t>(binom(i_midx[j], l_midx[j])) / static_cast<bry_float_t>(binom(to_degree, l_midx[j]));
        }
        return transformation_coeff;
    };
    return makeBigMatrix(to_degree, degree, makeCoeff);
}

template <std::size_t DIM>
Eigen::MatrixXd BRY::BernsteinBasis<DIM>::getInverseTransformationMatrix(bry_deg_t degree) {
    auto makeCoeff = [&] (const auto& i_midx, const auto& l_midx) -> bry_float_t {
        bry_float_t transformation_coeff = 1.0;

        for (bry_idx_t j = 0; j < DIM; ++j) {
            transformation_coeff *= static_cast<bry_float_t>(binom(degree - l_midx[j], degree - i_midx[j]));
        }

        // If the sum of all the i indices is even, the coefficient is positive, else negative
        bool neg = std::accumulate(i_midx.begin(), i_midx.end(), 0) % 2;
        return neg ? -transformation_coeff : transformation_coeff;
    };
    return makeBigMatrix(degree, degree, makeCoeff);
}

template <std::size_t DIM>
template <typename COEFF_LAM>
Eigen::MatrixXd BRY::BernsteinBasis<DIM>::makeBigMatrix(bry_deg_t to_degree, bry_deg_t from_degree, COEFF_LAM makeCoeff) {
    Eigen::MatrixXd matrix(pow(to_degree + 1, DIM), pow(from_degree + 1, DIM));
    matrix.setZero();

    // In place multi index for 'l' indices
    std::array<bry_idx_t, DIM> l_idx;
    std::array<bry_idx_t, DIM> i_idx;

    MultiIndex<ExhaustiveIncrementerWrap> i_midx(i_idx.data(), DIM, true, to_degree + 1);
    for (; !i_midx.last(); ++i_midx) {
        
        // Use the row midx as the bounds for the column (i) iterator
        std::vector<bry_idx_t> index_bounds;
        index_bounds.reserve(DIM);
        for (bry_idx_t r_i : i_midx)
            index_bounds.push_back(std::min(r_i + 1, from_degree + 1));

        // Create the column multi index
        MultiIndex<BoundedExhaustiveIncrementerWrap> l_midx(l_idx.data(), DIM, true, index_bounds, from_degree + 1);
        for (; !l_midx.last(); ++l_midx) {


            matrix(i_midx.inc().wrappedIdx(), l_midx.inc().wrappedIdx()) = makeCoeff(i_midx, l_midx);
        }
    }
    return matrix;
}

//Eigen::Tensor<bry_float_t, DIM> getDenominatorTensor()