#pragma once

#include "BernsteinTransform.h"
#include "MultiIndex.h"
#include "Operations.h"
#include "Types.h"

//template <std::size_t DIM>
//Polynomial<DIM, BRY::Basis::Bernstein> BRY::BernsteinBasisTransform<DIM>::to(const Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_int_t degree_increase = 0) {
//    
//}

template <std::size_t DIM>
BRY::Matrix BRY::BernsteinBasisTransform<DIM>::pwrToBernMatrix(bry_int_t degree, bry_int_t degree_increase) {
    bry_int_t to_degree = degree + degree_increase;
    auto makeCoeff = [&] (const auto& i_midx, const auto& l_midx) -> bry_float_t {
        // Compute the transformation coefficient in a numerically stable way
        bry_float_t transformation_coeff = 1.0;
        for (bry_int_t j = 0; j < DIM; ++j) {
            transformation_coeff *= static_cast<bry_float_t>(binom(i_midx[j], l_midx[j])) / static_cast<bry_float_t>(binom(to_degree, l_midx[j]));
        }
        return transformation_coeff;
    };
    return makeBigMatrix(to_degree, degree, makeCoeff);
}

template <std::size_t DIM>
BRY::Matrix BRY::BernsteinBasisTransform<DIM>::bernToPwrMatrix(bry_int_t degree) {
    auto makeCoeff = [&] (const auto& i_midx, const auto& l_midx) -> bry_float_t {
        bry_float_t transformation_coeff = 1.0;

        bool neg = false;
        for (bry_int_t j = 0; j < DIM; ++j) {
            transformation_coeff *= static_cast<bry_float_t>(binom(degree - l_midx[j], degree - i_midx[j]) * binom(degree, l_midx[j]));
            if ((i_midx[j] - l_midx[j]) % 2 != 0) 
                neg = !neg;
        }

        return neg ? -transformation_coeff : transformation_coeff;
    };
    return makeBigMatrix(degree, degree, makeCoeff);
}

template <std::size_t DIM>
template <typename COEFF_LAM>
BRY::Matrix BRY::BernsteinBasisTransform<DIM>::makeBigMatrix(bry_int_t to_degree, bry_int_t from_degree, COEFF_LAM makeCoeff) {
    Matrix matrix(pow(to_degree + 1, DIM), pow(from_degree + 1, DIM));
    matrix.setZero();

    // In place multi index for 'l' indices
    std::array<bry_int_t, DIM> l_idx;
    std::array<bry_int_t, DIM> i_idx;

    MultiIndex<ExhaustiveIncrementerWrap> i_midx(i_idx.data(), DIM, true, to_degree + 1);
    for (; !i_midx.last(); ++i_midx) {
        
        // Use the row midx as the bounds for the column (i) iterator
        std::vector<bry_int_t> index_bounds;
        index_bounds.reserve(DIM);
        for (bry_int_t r_i : i_midx)
            index_bounds.push_back(std::min(r_i + 1, from_degree + 1));

        // Create the column multi index
        MultiIndex<BoundedExhaustiveIncrementerWrap> l_midx(l_idx.data(), DIM, true, index_bounds, from_degree + 1);
        for (; !l_midx.last(); ++l_midx) {
            matrix(i_midx.inc().wrappedIdx(), l_midx.inc().wrappedIdx()) = makeCoeff(i_midx, l_midx);
        }
    }
    return matrix;
}