#pragma once

#include "Shrink.h"

#include "lemon/Logging.h"

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

template <std::size_t DIM>
static BRY::Polynomial<DIM, Basis::Power> BRY::Shrink::pruneHigherOrder(const Polynomial<DIM, Basis::Power>& p, bry_int_t new_degree) {
    std::array<bry_int_t, DIM> offsets = makeUniformArray<bry_int_t, DIM>(0);
    std::array<bry_int_t, 2> extents = makeUniformArray<bry_int_t, DIM>(new_degree);
    Eigen::Tensor<bry_float_t, DIM> shrunken_tensor = p.tensor().slice(offsets, extents);
    return BRY::Polynomial<DIM, Basis::Power>(std::move(shrunken_tensor));
}

template <std::size_t DIM>
static BRY::Polynomial<DIM, Basis::Power> BRY::Shrink::upperBound01(const Polynomial<DIM, Basis::Power>& p, bry_int_t new_degree) {
    
}

template <std::size_t DIM>
static BRY::Polynomial<DIM, Basis::Power> BRY::Shrink::lowerBound01(const Polynomial<DIM, Basis::Power>& p, bry_int_t new_degree) {

}