#pragma once

#include "BernsteinTransform.h"
#include "MultiIndex.h"
#include "Operations.h"
#include "Types.h"

#include "lemon/Logging.h"

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

template <std::size_t DIM>
std::pair<BRY::bry_float_t, bool> BRY::BernsteinBasisTransform<DIM>::infBound(const BRY::Polynomial<DIM, BRY::Basis::Bernstein>& p) {
    Eigen::Tensor<bry_float_t, 0> min = p.tensor().minimum();

    bry_float_t min_coeff = min();
    Eigen::Vector<bry_int_t, DIM> vertex_idx = Eigen::Vector<bry_int_t, DIM>::Zero();
    MultiIndex<ExhaustiveIncrementer> midx(vertex_idx.data(), DIM, true, 2);
    for (; !midx.last(); ++midx) {
        bry_float_t vertex_val = p.tensor()(p.degree() * vertex_idx);
        if (vertex_val == min_coeff) {
            return std::make_pair(min_coeff, true);
        }
    }

    return std::make_pair(min_coeff, false);
}

template <std::size_t DIM>
BRY::bry_float_t BRY::BernsteinBasisTransform<DIM>::infBoundGap(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, bool vertex_condition, bry_int_t degree_increase) {
    if (vertex_condition)
        return 0.0;

    bry_float_t epsilon = 0.0;
    std::array<bry_int_t, DIM> idx;
    MultiIndex<ExhaustiveIncrementer> midx(idx.data(), DIM, true, p.degree() + 1);
    for (; !midx.last(); ++midx) {
        bry_int_t multiplier = 0;
        for (bry_int_t d = 0; d < DIM; ++d) {
            if (midx[d] != 0) {
                multiplier += (midx[d] - 1) * (midx[d] - 1);
            }
        }
        //DEBUG("midx: " << midx << " adding " << static_cast<bry_float_t>(multiplier) * std::abs(p.coeff(idx)) << " from coeff: " << (p.coeff(idx)));
        epsilon += static_cast<bry_float_t>(multiplier) * std::abs(p.coeff(idx));
    }
    //DEBUG("epsilon: " << epsilon << " degree: " << p.degree());
    //PAUSE;

    bry_float_t raised_deg_f = static_cast<bry_float_t>(p.degree() + degree_increase);
    return epsilon * (raised_deg_f - 1.0) / (raised_deg_f * raised_deg_f);
}