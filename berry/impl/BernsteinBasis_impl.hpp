#pragma once

#include "Bernstein.h"
#include "MultiIndex.h"

template <std::size_t DIM>
Polynomial<DIM, BRY::Basis::Bernstein> BRY::BernsteinBasis<DIM>::to(const Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_deg_t degree_increase = 0) {
    
}

template <std::size_t DIM>
Eigen::Tensor<bry_float_t, 2> getTransformationMatrix(const Polynomial<DIM, Basis::Power>& p, bry_deg_t degree_increase) {
    Eigen::Tensor<bry_float_t, 2> matrix(p.nMonomails(), p.nMonomails());
    matrix.setZero();
    for (auto row_midx = mIdx(DIM, p.degree() + 1); !row_midx.last(); ++row_midx) {
    }
}

Eigen::Tensor<bry_float_t, DIM> getDenominatorTensor()