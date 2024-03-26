#pragma once

#include <cmath>
//#include <experimental/array>

#include "Polynomial.h"

template <std::size_t DIM>
BRY::Polynomial<DIM>::Polynomial(bry_deg_t degree)
    : m_degree(degree), m_container(std::pow(degree + 1, DIM), bry_float_t{})
{ }

template <std::size_t DIM>
std::size_t BRY::Polynomial<DIM>::degree() const {
    return m_degree;
}

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
//template <std::size_t DIM>
//template <typename ... FLTS>
//BRY::bry_float_t BRY::Polynomial<DIM>::operator()(FLTS ... x) const {
//    static_assert(is_uniform_convertible_type<bry_float_t, FLTS ...>(), "All parameters passed to `operator()` must be float type (`bry_float_t`)");
//    static_assert(sizeof...(FLTS) == DIM, "Number of x parameters must match the dimension of the polynomial");
//    auto x_arr = makeArray<bry_float_t>(x ...);
//
//    //for (bry_float_t coeff)
//}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM>& p) {
    bool first = true;
    for (std::size_t i = 0; i < p.m_container.size(); ++i) {
        if (p.m_container[i] == BRY::bry_float_t{}) 
            continue;

        if (!first)
            os << BRY_LOG_WHITE(" + ");
        first = false;

        os << BRY_LOG_BYELLOW(p.m_container[i]);
        std::array<BRY::bry_deg_t, DIM> degrees = p.unwrap(i);
        for (std::size_t dim = 0; dim < DIM; ++dim) {
            if (degrees[dim] > 0)
                os << BRY_LOG_WHITE("(x" << dim << "^") << BRY_LOG_BGREEN(degrees[dim]) << BRY_LOG_WHITE(")");
        }
    }
    return os;
}

template <std::size_t DIM>
std::size_t BRY::Polynomial<DIM>::wrap(const std::array<bry_deg_t, DIM>& exponents) const {
    std::size_t ret_idx = 0;

    for (std::size_t i = 0; i < DIM; ++i) {

#ifdef BRY_ENABLE_BOUNDS_CHECK
        ASSERT(exponents[i] <= m_degree, "Degree of x" << i << " (" << exponents[i] << ") is out of bounds (degree of polynomial is " << m_degree << ")");
#endif

        ret_idx += std::pow(m_degree + 1, i) * exponents[i];
    }
    return ret_idx;
    
}

template <std::size_t DIM>
std::array<BRY::bry_deg_t, DIM> BRY::Polynomial<DIM>::unwrap(std::size_t idx) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(idx < std::pow(m_degree + 1, DIM), "Wrapped index (" << idx << ") is out of bounds (max index is " << std::pow(m_degree + 1, DIM) << ")");
#endif

    std::array<bry_deg_t, DIM> exponents;

    std::size_t temp_mod = idx;
    for (std::size_t i = 0; i < DIM; ++i) {
        std::size_t temp = std::pow(m_degree + 1, DIM - i - 1);
        exponents[DIM - i - 1] = (temp_mod - idx % temp) / temp;
        temp_mod = idx % temp;
    }
    return exponents;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator+(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2) {
    const BRY::Polynomial<DIM>* p_big;
    const BRY::Polynomial<DIM>* p_small;
    if ((p_1.m_degree > p_2.m_degree)) {
        p_big = &p_1;
        p_small = &p_2;
    } else {
        p_big = &p_2;
        p_small = &p_1;
    }

    BRY::Polynomial<DIM> p_new = *p_big;

    std::size_t small_degree_powers[DIM];
    std::size_t big_degree_powers[DIM];
    for (std::size_t d = 0; d < DIM; ++d) {
        small_degree_powers[d] = std::pow(p_small->m_degree, d);
        big_degree_powers[d] = std::pow(p_big->m_degree, d);
    }

    std::size_t deg_diff = p_big->m_degree - p_small->m_degree;
    std::size_t i_big = 0;

    for (std::size_t i_small = 0; i_small < p_small->m_container.size(); ++i_small) {
        //DEBUG("B4: i_small: " << i_small << ", i_big: " << i_big);
        p_new.m_container[i_big] += p_small->m_container[i_small];
        
        i_big += 1;
        for (int64_t d = DIM - 1; d >= 0; --d) {
            //DEBUG("d: " << d);
            if (((i_small + 1) % small_degree_powers[d]) == 0) {
                i_big += big_degree_powers[d] * deg_diff;
            }
        }
        //DEBUG("AF: i_small: " << i_small << ", i_big: " << i_big);
    }

    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator-(const BRY::Polynomial<DIM>& p) {
    return -1 * p;
}


template <std::size_t DIM>
BRY::Polynomial<DIM> operator-(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2) {
    return p_1 + -1 * p_2;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(BRY::bry_float_t scalar, const BRY::Polynomial<DIM>& p) {
    BRY::Polynomial<DIM> p_new = p;
    for (BRY::bry_float_t& c : p_new.m_container)
        c = scalar * c;
    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(const BRY::Polynomial<DIM>& p, BRY::bry_float_t scalar) {
    return scalar * p;
}

/* TODO: Make this faster */
template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2) {

    const BRY::Polynomial<DIM>* p_big;
    const BRY::Polynomial<DIM>* p_small;
    if ((p_1.m_degree > p_2.m_degree)) {
        p_big = &p_1;
        p_small = &p_2;
    } else {
        p_big = &p_2;
        p_small = &p_1;
    }

    BRY::Polynomial<DIM> p_new(p_1.m_degree + p_2.m_degree);

    for (std::size_t i = 0; i < p_1.m_container.size(); ++i) {
        std::array<BRY::bry_deg_t, DIM> i_unwp = p_1.unwrap(i);
        for (std::size_t j = 0; j < p_2.m_container.size(); ++j) {
            std::array<BRY::bry_deg_t, DIM> j_unwp = p_2.unwrap(j);
            std::array<BRY::bry_deg_t, DIM> k_unwp;
            for (std::size_t d = 0; d < DIM; ++d)
                k_unwp[d] = i_unwp[d] + j_unwp[d];
            p_new.m_container[p_new.wrap(k_unwp)] += p_1.m_container[i] + p_2.m_container[j];
        }
    }

    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator^(const BRY::Polynomial<DIM>& p_1, BRY::bry_deg_t deg) {

}