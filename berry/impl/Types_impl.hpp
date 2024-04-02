#pragma once

#include "Types.h"

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
Eigen::Tensor<BRY::bry_float_t, DIM> BRY::makeIncrementTensor(const std::array<bry_deg_t, DIM>& dimensions, bry_deg_t increment_idx, bry_deg_t offset) {
    std::array<bry_deg_t, DIM> before_broadcast_dims = makeUniformArray<bry_deg_t, DIM>(1);
    before_broadcast_dims[increment_idx] = dimensions[increment_idx];
    Eigen::Tensor<bry_float_t, DIM> increment_vector(before_broadcast_dims);

    std::array<bry_deg_t, DIM> midx = makeUniformArray<bry_deg_t, DIM>(0);
    for (bry_deg_t i = 0; i < dimensions[increment_idx]; ++i) {
        ++midx[increment_idx] = i;
        increment_vector(midx) = static_cast<bry_float_t>(i + offset);
    }

    std::array<bry_deg_t, DIM> bcast = dimensions;
    bcast[increment_idx] = 1;

    return increment_vector.broadcast(bcast);
}