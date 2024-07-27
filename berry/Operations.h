#pragma once

#include "Options.h"
#include "Types.h"
#include "MultiIndex.h"

#include <vector>

namespace BRY {

static BRY_INL std::size_t factorial(std::size_t n);

static BRY_INL std::size_t binom(std::size_t n, std::size_t k);

/// @brief Creates a vector of all binomial coefficients (n'th row of Pascal's triangle)
/// of the form [(n choose 0), (n choose 1), ..., (n choose n)]
/// @param n 
/// @return Array of size `n + 1` of all binomial coefficients from k = 0 to k = n
static std::vector<std::size_t> pascalRow(std::size_t n);

template <class ITERABLE_T>
static BRY_INL std::size_t multinom(std::size_t n, const ITERABLE_T& multi_k);

static BRY_INL std::size_t multinom(const MultiIndex<>& idx);

/// @brief Exponent for integers and (smallish) positive exponents
/// @param x Base
/// @param n Exponent (positive)
/// @return Integer x^n
static BRY_INL bry_int_t pow(bry_int_t x, bry_int_t n);

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

template <std::size_t DIM>
static Matrix makeDegreeChangeTransform(bry_int_t from_deg, bry_int_t to_deg);

}

#include "impl/Operations_impl.hpp"