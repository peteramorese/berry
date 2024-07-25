#pragma once

#include "Options.h"
#include "Types.h"
#include "Polynomial.h"

#include <vector>
#include <array>
#include <memory>

namespace BRY {

namespace Shrink {

/// @brief Reduce degree of polynomial by simply pruning higher order terms
/// @tparam DIM Dimension of polynomial
/// @param p Polynomial to shrink
/// @param new_degree Degree of the shrunken polynomial
/// @return Reduced degree polynomial
template <std::size_t DIM>
static Polynomial<DIM, Basis::Power> pruneHigherOrder(const Polynomial<DIM, Basis::Power>& p, bry_int_t new_degree);

/// @brief Reduce degree of polynomial such that the lower order polynomial upper bounds the original on the hypercube `[0, 1]^DIM`
/// @tparam DIM Dimension of polynomial
/// @param p Polynomial to shrink
/// @param new_degree Degree of the shrunken polynomial
/// @return Reduced degree polynomial (upper bound of `p`)
template <std::size_t DIM>
static Polynomial<DIM, Basis::Power> upperBound01(const Polynomial<DIM, Basis::Power>& p, bry_int_t new_degree);

/// @brief Reduce degree of polynomial such that the lower order polynomial lower bounds the original on the hypercube `[0, 1]^DIM`
/// @tparam DIM Dimension of polynomial
/// @param p Polynomial to shrink
/// @param new_degree Degree of the shrunken polynomial
/// @return Reduced degree polynomial (lower bound of `p`)
template <std::size_t DIM>
static Polynomial<DIM, Basis::Power> lowerBound01(const Polynomial<DIM, Basis::Power>& p, bry_int_t new_degree);

}

}

#include "impl/Shrink_impl.hpp"