#pragma once

#include <cmath>
//#include <experimental/array>

#include "Polynomial.h"

template <std::size_t DIM>
BRY::Polynomial<DIM>::Polynomial(bry_deg_t degree)
    : m_degree(degree), m_container(std::pow(degree, DIM), bry_float_t{})
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

/* TODO */
template <std::size_t DIM>
template <typename ... FLTS>
BRY::bry_float_t BRY::Polynomial<DIM>::operator()(FLTS ... x) const {
    static_assert(is_uniform_convertible_type<bry_float_t, FLTS ...>(), "All parameters passed to `operator()` must be float type (`bry_float_t`)");
    static_assert(sizeof...(FLTS) == DIM, "Number of x parameters must match the dimension of the polynomial");
    auto x_arr = makeArray<bry_float_t>(x ...);

    //for (bry_float_t coeff)
}

template <std::size_t DIM>
std::ostream& BRY::Polynomial<DIM>::print(std::ostream& os) const {
    bool first = true;
    for (std::size_t i = 0; i < m_container.size(); ++i) {
        if (m_container[i] == BRY::bry_float_t{}) 
            continue;

        if (!first)
            os << BRY_LOG_WHITE(" + ");
        first = false;

        os << BRY_LOG_BWHITE(m_container[i]);
        std::array<BRY::bry_deg_t, DIM> degrees = unwrap(i);
        for (std::size_t dim = 0; dim < DIM; ++dim) {
            if (degrees[dim] > 0)
                os << BRY_LOG_GREEN("x") << BRY_LOG_BGREEN(degrees[dim]);
        }
    }
    return os;
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

template <std::size_t DIM>
BRY::Polynomial<DIM> BRY::operator+(const Polynomial<DIM>& p_1, const Polynomial<DIM>& p_2) {
    const Polynomial<DIM>* p_big, p_small;
    if ((p_1.m_degree > p_2.m_degree)) {
        p_big = &p_1;
        p_small = &p_2;
    } else {
        p_big = &p_2;
        p_small = &p_1;
    }

    Polynomial<DIM> p_new = *p_big;

    for (std::size_t i = 0; i < p_small->container.size(); ++i) {
        p_new.m_container[i] += p_small->conatiner[i];
    }
    
    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> BRY::operator-(const Polynomial<DIM>& p_1, const Polynomial<DIM>& p_2) {

}

template <std::size_t DIM>
BRY::Polynomial<DIM> BRY::operator*(const Polynomial<DIM>& p_1, const Polynomial<DIM>& p_2) {

}

template <std::size_t DIM>
BRY::Polynomial<DIM> BRY::operator^(const Polynomial<DIM>& p_1, bry_deg_t deg) {

}