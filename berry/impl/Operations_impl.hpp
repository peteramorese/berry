#pragma once

#include "Operations.h"

#include <cmath>
#include <algorithm>
#include <numeric>

std::size_t BRY::factorial(std::size_t n) {
    std::size_t val(1.0);
    for (std::size_t i = 2; i < n; ++i)
        val *= i;
    return val;
}

std::size_t BRY::binom(std::size_t n, std::size_t k) {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(k <= n, "`k` must be <= `n`");
#endif

    std::size_t val(1);
    for (std::size_t i = 1; i < k + 1; ++i) {
        val *= (n + 1 - i);
        val /= i;
    }
    return val;
}

std::vector<std::size_t> BRY::pascalRow(std::size_t n) {
    std::vector<std::size_t> vals(n + 1);
    vals.front() = 1;
    for (std::size_t i = 1; i < n + 1; ++i) {
        vals[i] = (vals[i - 1] * (n + 1 - i)) / i;
    }
    return vals;
}

template <class ITERABLE_T>
std::size_t BRY::multinom(std::size_t n, const ITERABLE_T& multi_k) {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(std::accumulate(multi_k.begin(), multi_k.end(), 0) == n, "the sum of `k` must be == `n`");
#endif

    std::size_t val(1.0);
    auto it = multi_k.begin();
    std::size_t k_partial_sum = *(it++);
    for (; it != multi_k.end(); ++it) {
        k_partial_sum += *it;
        val *= binom(k_partial_sum, *it);
    }
    return val;
}

std::size_t BRY::multinom(const MultiIndex<>& idx) {
    return multinom(idx.inc().indexConstraint(), idx);
}

BRY::bry_int_t BRY::pow(bry_int_t x, std::size_t n) {
    if (n == 0)
        return 1;

    for (std::size_t i = 1; i < n; ++i)
        x *= x;
    return x;
}