#pragma once

#include <cmath>
//#include <experimental/array>

#include "Polynomial.h"

template <std::size_t DIM>
BRY::Polynomial<DIM>::Polynomial(bry_deg_t degree)
    : m_container(degree, bry_float_t{})
{ }

template <std::size_t DIM>
template <typename ... DEGS>
BRY::bry_float_t& BRY::Polynomial<DIM>::coeff(DEGS ... exponents) {
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_container[wrap(makeArray<bry_deg_t>(exponents...))];
}

template <std::size_t DIM>
template <typename ... DEGS>
const BRY::bry_float_t& BRY::Polynomial<DIM>::coeff(DEGS ... exponents) const {
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_container[wrap(makeArray<bry_deg_t>(exponents...))];
}

template <std::size_t DIM>
std::size_t BRY::Polynomial<DIM>::wrap(const std::array<bry_deg_t, DIM>& exponents) const {
    std::size_t ret_idx = 0;

    for (std::size_t i = 0; i < DIM; ++i) {

#ifdef BRY_ENABLE_BOUNDS_CHECK
        ASSERT(exponents[i] < m_degree, "Degree of x" << i << " (" << exponents[i] << ") is out of bounds (degree of polynomial is " << m_degree << ")");
#endif

        ret_idx += std::pow(m_degree, i) * exponents[i];
    }
    return ret_idx;
    
}

template <std::size_t DIM>
std::array<BRY::bry_deg_t, DIM> BRY::Polynomial<DIM>::unwrap(std::size_t idx) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(idx < std::pow(m_degree, DIM), "Wrapped index (" << idx << ") is out of bounds (max index is " << std::pow(m_degree, DIM) << ")");
#endif

    std::array<bry_deg_t, DIM> exponents;

    std::size_t temp_mod = idx;
    for (std::size_t i = 0; i < DIM; ++i) {
        std::size_t temp = std::pow(m_degree, DIM - i - 1);
        exponents[DIM - i - 1] = (temp_mod - idx % temp) / temp;
        temp_mod = idx % temp;
    }
    return exponents;
}