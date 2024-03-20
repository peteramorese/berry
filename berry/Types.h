#pragma once

#include <type_traits>
#include <array>
#include <tuple>

namespace BRY {

/* Type aliases */

/// @brief General floating point
typedef double bry_float_t;

/// @brief General index
typedef uint32_t bry_idx_t;

/// @brief Degree of polynomials
typedef uint64_t bry_deg_t;


/* Helpful functions*/

template <typename T, typename... ARGS_T>
constexpr bool is_uniform_type();

template <typename T, typename... ARGS_T>
constexpr bool is_uniform_convertible_type();

template <typename T, typename... ARGS_T>
std::array<T, sizeof...(ARGS_T)> makeArray(ARGS_T&&... args);

}

#include "impl/Types_impl.hpp"
