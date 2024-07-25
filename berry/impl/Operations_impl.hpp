#pragma once

#include "Operations.h"
#include "MultiIndex.h"

#include "lemon/Logging.h"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace _BRY {

template <std::size_t I, typename T, typename ... ARGS_T>
constexpr bool _is_uniform_type() {
    using T_I = std::tuple_element<I, std::tuple<ARGS_T...>>::type;
    if constexpr (std::is_convertible<T_I, T>::value) {
        if constexpr (I < sizeof...(ARGS_T) - 1) {
            return _is_uniform_type<I + 1, T, ARGS_T...>();
        } else {
            return true;
        }
    } else {
        return false;
    }
}

template <std::size_t I, typename T, typename... ARGS_T>
constexpr bool _is_uniform_convertible_type() {
    using T_I = std::tuple_element<I, std::tuple<ARGS_T...>>::type;
    if constexpr (std::is_convertible<T_I, T>::value) {
        if constexpr (I < sizeof...(ARGS_T) - 1) {
            return _is_uniform_convertible_type<I + 1, T, ARGS_T...>();
        } else {
            return true;
        }
    } else {
        return false;
    }
}

template <std::size_t I, typename T, typename... ARGS_T>
void _setArrayElement(std::array<T, sizeof...(ARGS_T)>& array, const std::tuple<ARGS_T&&...>& args_tuple) {
    if constexpr (I < sizeof...(ARGS_T)) {
        std::get<I>(array) = std::get<I>(args_tuple);
        _setArrayElement<I + 1, T, ARGS_T...>(array, args_tuple);
    }
}

template <std::size_t I, typename... ARGS_T>
void _setExponentVecElement(BRY::ExponentVec<sizeof...(ARGS_T)>& exp, const std::tuple<ARGS_T&&...>& args_tuple) {
    if constexpr (I < sizeof...(ARGS_T)) {
        exp[I] = std::get<I>(args_tuple);
        _setExponentVecElement<I + 1, ARGS_T...>(exp, args_tuple);
    }
}

}

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

template <typename T, typename... ARGS_T>
constexpr bool BRY::is_uniform_type() {
    return _BRY::_is_uniform_type<0, T, ARGS_T...>();
}

template <typename T, typename... ARGS_T>
constexpr bool BRY::is_uniform_convertible_type() {
    return _BRY::_is_uniform_convertible_type<0, T, ARGS_T...>();
}

template <typename T, typename... ARGS_T>
std::array<T, sizeof...(ARGS_T)> BRY::makeArray(ARGS_T&&... args) {
    static_assert(is_uniform_convertible_type<T, ARGS_T ...>(), "All parameters passed to `makeArray` must be of the same type");
    std::tuple<ARGS_T&&...> args_tuple(std::forward<ARGS_T>(args)...);

    std::array<T, sizeof...(ARGS_T)> array;
    _BRY::_setArrayElement<0, T, ARGS_T...>(array, args_tuple);

    return array;
}

template <typename T, std::size_t SZ>
std::array<T, SZ> BRY::makeUniformArray(const T& fill_val) {
    std::array<T, SZ> array;
    array.fill(fill_val);
    return array;
}

template <typename... ARGS_T>
BRY::ExponentVec<sizeof...(ARGS_T)> BRY::makeExponentVec(ARGS_T&&... args) {
    static_assert(is_uniform_convertible_type<bry_float_t, ARGS_T ...>(), "All parameters passed to `makeExponentVec` must be `bry_float_t` type");
    std::tuple<ARGS_T&&...> args_tuple(std::forward<ARGS_T>(args)...);

    ExponentVec<sizeof...(ARGS_T)> exp;
    _BRY::_setExponentVecElement<0, ARGS_T...>(exp, args_tuple);

    return exp;
}

template <std::size_t DIM>
Eigen::Tensor<BRY::bry_float_t, DIM> BRY::makeIncrementTensor(const std::array<bry_int_t, DIM>& dimensions, bry_int_t increment_idx, bry_int_t offset) {
    std::array<bry_int_t, DIM> before_broadcast_dims = makeUniformArray<bry_int_t, DIM>(1);
    before_broadcast_dims[increment_idx] = dimensions[increment_idx];
    Eigen::Tensor<bry_float_t, DIM> increment_vector(before_broadcast_dims);

    std::array<bry_int_t, DIM> midx = makeUniformArray<bry_int_t, DIM>(0);
    for (bry_int_t i = 0; i < dimensions[increment_idx]; ++i) {
        ++midx[increment_idx] = i;
        increment_vector(midx) = static_cast<bry_float_t>(i + offset);
    }

    std::array<bry_int_t, DIM> bcast = dimensions;
    bcast[increment_idx] = 1;

    return increment_vector.broadcast(bcast);
}

template <std::size_t DIM>
static BRY::Matrix BRY::makeDegreeChangeTransform(bry_int_t from_deg, bry_int_t to_deg) {
    bry_int_t rows = pow(to_deg + 1, DIM);
    bry_int_t cols = pow(from_deg + 1, DIM);
    //DEBUG("rows: " << rows << " cols: " << cols);
    Matrix tf = Matrix::Zero(rows, cols);
    auto row_midx = mIdxW(DIM, to_deg + 1);
    auto col_midx = mIdxW(DIM, from_deg + 1);
    while (!row_midx.last() && !col_midx.last()) {
        bry_int_t i = 0;
        //DEBUG("row: " << row_midx << " col: " << col_midx);
        for (; i < DIM; ++i) {
            if (row_midx[i] < col_midx[i]) {
                ++col_midx;
                break;
            } else if (row_midx[i] > col_midx[i]) {
                ++row_midx;
                break;
            }
        }
        if (i < DIM - 1)
            continue;
        //DEBUG("accessing " << row_midx.inc().wrappedIdx() << ", " << col_midx.inc().wrappedIdx());
        tf(row_midx.inc().wrappedIdx(), col_midx.inc().wrappedIdx()) = 1.0;
        ++row_midx;
        ++col_midx;
    }
    return tf;
}