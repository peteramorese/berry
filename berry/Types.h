#pragma once

#include <type_traits>
#include <array>
#include <tuple>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

namespace BRY {

/* Type aliases */

/// @brief General floating point
typedef double bry_float_t;

/// @brief General complex floating point
typedef std::complex<bry_float_t> bry_complex_t;

/// @brief General integer
typedef int64_t bry_int_t;

/// @brief Vector of multinomial exponents
template <std::size_t DIM>
using ExponentVec = Eigen::Vector<bry_int_t, DIM>;

using Vector = Eigen::Matrix<bry_float_t, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<bry_float_t, Eigen::Dynamic, Eigen::Dynamic>;

}
