#pragma once

#include "Types.h"

#include <tuple>

namespace _BRY {

template <uint64_t I, typename T, typename ... ARGS_T>
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

template <uint64_t I, typename T, typename ... ARGS_T>
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

}

template <typename T, typename ... ARGS_T>
constexpr bool BRY::is_uniform_type() {
    return _BRY::_is_uniform_type<0, T, ARGS_T...>();
}

template <typename T, typename ... ARGS_T>
constexpr bool BRY::is_uniform_convertible_type() {
    return _BRY::_is_uniform_convertible_type<0, T, ARGS_T...>();
}
