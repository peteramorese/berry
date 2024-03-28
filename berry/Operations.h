#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"
#include "MultiIndex.h"

#include <vector>

namespace BRY {

static BRY_INL BRY::bry_float_t factorial(std::size_t n);

static BRY_INL BRY::bry_float_t binom(std::size_t n, std::size_t k);

template <class ITERABLE_T>
static BRY_INL BRY::bry_float_t multinom(std::size_t n, const ITERABLE_T& multi_k);

static BRY_INL BRY::bry_float_t multinom(const MultiIndex& idx);

}

#include "impl/Operations_impl.hpp"