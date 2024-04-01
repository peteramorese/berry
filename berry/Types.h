#pragma once

#include <type_traits>
#include <array>
#include <tuple>

#include <Eigen/Core>

namespace BRY {

/* Type aliases */

/// @brief General floating point
typedef double bry_float_t;

/// @brief General floating point
typedef std::complex<bry_float_t> bry_complex_t;

/// @brief General index
typedef int64_t bry_idx_t;

/// @brief Degree of polynomials
typedef int64_t bry_deg_t;

/// @brief Vector of multinomial exponents
template <std::size_t DIM>
using ExponentVec = Eigen::Vector<bry_deg_t, DIM>;


/* Helpful functions*/

template <typename T, typename... ARGS_T>
constexpr bool is_uniform_type();

template <typename T, typename... ARGS_T>
constexpr bool is_uniform_convertible_type();

template <typename T, typename... ARGS_T>
std::array<T, sizeof...(ARGS_T)> makeArray(ARGS_T&&... args);

template <typename T, std::size_t SZ>
std::array<T, SZ> makeUniformArray(const T& fill_val);

template <typename... ARGS_T>
ExponentVec<sizeof...(ARGS_T)> makeExponentVec(ARGS_T&&... args);

}

#include "impl/Types_impl.hpp"
