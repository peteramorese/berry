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