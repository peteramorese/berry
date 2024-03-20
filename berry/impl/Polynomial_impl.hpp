#pragma once

#include <cmath>

#include "Polynomial.h"

template <std::size_t DIM, class DEGREE_T>
BRY::Polynomial<DIM, DEGREE_T>::Polynomial(std::size_t degree)
    : m_degree(degree)
{}

template <std::size_t DIM, class DEGREE_T>
template <typename ... DEGS>
BRY::bry_float_t BRY::Polynomial<DIM, DEGREE_T>::coeff(DEGS ... degrees) {
    static_assert(is_uniform_convertible_type<bry_idx_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of degrees must match the dimension of the polynomial");

}

template <std::size_t DIM, class DEGREE_T>
std::size_t BRY::Polynomial<DIM, DEGREE_T>::wrap(const std::array<bry_deg_t, DIM>& degrees) const {
    std::size_t ret_idx = 0;

    for (std::size_t i = 0; i < DIM; ++i) {

#ifdef BRY_ENABLE_BOUNDS_CHECK
        ASSERT(degrees[i] < m_degree, "Degree of x" << i << " (" << degrees[i] << ") is out of bounds (degree of polynomial is " << m_degree << ")");
#endif

        ret_idx += std::pow(m_degree, i) * degrees[i];
    }
    return ret_idx;
    
}

template <std::size_t DIM, class DEGREE_T>
std::array<BRY::bry_deg_t, DIM> BRY::Polynomial<DIM, DEGREE_T>::unwrap(std::size_t idx) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    //ASSERT(idx < m_container.arr.size(), "Wrapped index (" << idx << ") is out of bounds (max index is " << m_container.arr.size() << ")");
#endif

    std::array<bry_deg_t, DIM> degrees;

    std::size_t temp_mod = idx;
    for (std::size_t i = 0; i < DIM; ++i) {
        std::size_t temp = std::pow(m_degree, DIM - i - 1);
        degrees[DIM - i - 1] = (temp_mod - idx % temp) / temp;
        temp_mod = idx % temp;
    }
    return degrees;
}