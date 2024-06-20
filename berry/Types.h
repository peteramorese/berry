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


/* Helpful functions*/

template <typename T, typename... ARGS_T>
constexpr bool is_uniform_type();

template <typename T, typename... ARGS_T>
constexpr bool is_uniform_convertible_type();

template <typename T, typename... ARGS_T>
static std::array<T, sizeof...(ARGS_T)> makeArray(ARGS_T&&... args);

template <typename T, std::size_t SZ>
static std::array<T, SZ> makeUniformArray(const T& fill_val);

template <typename... ARGS_T>
static ExponentVec<sizeof...(ARGS_T)> makeExponentVec(ARGS_T&&... args);

template <std::size_t DIM>
static BRY_INL Eigen::Tensor<bry_float_t, DIM> makeIncrementTensor(const std::array<bry_int_t, DIM>& dimensions, bry_int_t increment_idx, bry_int_t offset = 0);

}

#include "impl/Types_impl.hpp"
