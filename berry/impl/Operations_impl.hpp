#pragma once

#include "Operations.h"

#include <cmath>
#include <algorithm>
#include <numeric>

BRY::bry_float_t BRY::factorial(std::size_t n) {
    BRY::bry_float_t val(1.0);
    for (std::size_t i = 2; i < n; ++i)
        val *= i;
    return val;
}

BRY::bry_float_t BRY::binom(std::size_t n, std::size_t k) {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(k <= n, "`k` must be <= `n`");
#endif

    BRY::bry_float_t val(1.0);
    for (std::size_t i = 1; i < k + 1; ++i) {
        val *= static_cast<BRY::bry_float_t>(n + 1 - i) / static_cast<BRY::bry_float_t>(i);
    }
    return val;
}

template <class ITERABLE_T>
BRY::bry_float_t BRY::multinom(std::size_t n, const ITERABLE_T& multi_k) {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(std::accumulate(multi_k.begin(), multi_k.end(), 0) == n, "the sum of `k` must be == `n`");
#endif

    BRY::bry_float_t val(1.0);
    auto it = multi_k.begin();
    std::size_t k_partial_sum = *(it++);
    for (; it != multi_k.end(); ++it) {
        k_partial_sum += *it;
        val *= binom(k_partial_sum, *it);
    }
    return val;
}

BRY::bry_float_t BRY::multinom(const MultiIndex& idx) {
    return multinom(idx.l1Norm(), idx);
}