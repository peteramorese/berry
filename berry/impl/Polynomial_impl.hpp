#pragma once

#include <cmath>
//#include <experimental/array>

#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"

template <std::size_t DIM>
BRY::Polynomial<DIM>::Polynomial(bry_deg_t degree)
    : m_tensor(makeUniformArray<bry_deg_t, DIM>(degree + 1))
{
    m_tensor.setZero();
}

template <std::size_t DIM>
BRY::Polynomial<DIM>::Polynomial(const Eigen::Tensor<bry_float_t, DIM>& tensor) 
    : m_tensor(tensor)
{}

template <std::size_t DIM>
BRY::Polynomial<DIM>::Polynomial(Eigen::Tensor<bry_float_t, DIM>&& tensor) 
    : m_tensor(std::move(tensor))
{}

template <std::size_t DIM>
std::size_t BRY::Polynomial<DIM>::degree() const {
    return m_tensor.dimension(0) - 1;
}

template <std::size_t DIM>
template <typename ... DEGS>
BRY::bry_float_t& BRY::Polynomial<DIM>::coeff(DEGS ... exponents) {
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_tensor(exponents...);
}

template <std::size_t DIM>
template <typename ... DEGS>
const BRY::bry_float_t& BRY::Polynomial<DIM>::coeff(DEGS ... exponents) const {
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_tensor(exponents...);
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
    std::array<BRY::bry_deg_t, DIM> idx_arr = BRY::makeUniformArray<BRY::bry_deg_t, DIM>(BRY::bry_deg_t{});
    std::size_t d = 0;
    bool first = true;

    auto iterate = [&] () {
        ++idx_arr[0];
        for (std::size_t i = 0; i < DIM - 1; ++ i) {
            if (idx_arr[i] > p.degree()) {
                idx_arr[i] = 0;
                ++idx_arr[i + 1];
            }
        }
    };

    while (idx_arr.back() <= p.degree()) {

        BRY::bry_float_t coeff = p.m_tensor(idx_arr);
        if (coeff == BRY::bry_deg_t{}) {
            iterate();
            continue;
        }

        if (!first)
            os << BRY_LOG_WHITE(" + ");
        first = false;

        os << BRY_LOG_BYELLOW(coeff);
        for (std::size_t dim = 0; dim < DIM; ++dim) {
            if (idx_arr[dim] > 0)
                os << BRY_LOG_WHITE("(x" << dim << "^") << BRY_LOG_BGREEN(idx_arr[dim]) << BRY_LOG_WHITE(")");
        }
        iterate();
    }
    return os;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator+(BRY::bry_float_t scalar, const BRY::Polynomial<DIM>& p) {
    BRY::Polynomial<DIM> p_new = p;
    *p_new.m_tensor.data() += scalar;
    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator+(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2) {
    const BRY::Polynomial<DIM>* p_big;
    const BRY::Polynomial<DIM>* p_small;
    if ((p_1.degree() > p_2.degree())) {
        p_big = &p_1;
        p_small = &p_2;
    } else {
        p_big = &p_2;
        p_small = &p_1;
    }

    std::array<std::pair<BRY::bry_idx_t, BRY::bry_idx_t>, DIM> paddings;
    for (std::pair<BRY::bry_idx_t, BRY::bry_idx_t>& pads : paddings) {
        pads.first = 0;
        pads.second = p_big->degree() - p_small->degree();
    }

    Eigen::Tensor<BRY::bry_float_t, DIM> new_tensor = p_small->m_tensor.pad(paddings);
    new_tensor += p_big->m_tensor;
    BRY::Polynomial<DIM> p_new(std::move(new_tensor));
    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator-(const BRY::Polynomial<DIM>& p) {
    return BRY::Polynomial<DIM>(-p.m_tensor);
}


template <std::size_t DIM>
BRY::Polynomial<DIM> operator-(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2) {
    return p_1 + -p_2;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(BRY::bry_float_t scalar, const BRY::Polynomial<DIM>& p) {
    return BRY::Polynomial<DIM>(scalar * p.m_tensor);
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(const BRY::Polynomial<DIM>& p, BRY::bry_float_t scalar) {
    return scalar * p;
}

/* TODO: Make this faster */
template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2) {

    //BRY::Polynomial<DIM> p_new(p_1.m_degree + p_2.m_degree);

    //for (std::size_t i = 0; i < p_1.m_container.size(); ++i) {
    //    BRY::ExponentVec<DIM> i_unwp = p_1.unwrap(i);
    //    for (std::size_t j = 0; j < p_2.m_container.size(); ++j) {
    //        BRY::ExponentVec<DIM> j_unwp = p_2.unwrap(j);
    //        BRY::ExponentVec<DIM> k_unwp;
    //        for (std::size_t d = 0; d < DIM; ++d)
    //            k_unwp[d] = i_unwp[d] + j_unwp[d];
    //        p_new.m_container[p_new.wrap(k_unwp)] += p_1.m_container[i] * p_2.m_container[j];
    //    }
    //}

    //return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator^(const BRY::Polynomial<DIM>& p, BRY::bry_deg_t deg) {
    BRY::Polynomial<DIM> p_new(p.m_degree * deg);

    for (BRY::MultiIndex midx(p.m_container.size(), deg); !midx.right(); ++midx) {
        BRY::ExponentVec<DIM> sum_exp = BRY::ExponentVec<DIM>::Zero();

        BRY::bry_float_t coeff(1.0);
        for (std::size_t monom_i = 0; monom_i < midx.size(); ++monom_i) {
            sum_exp += midx[monom_i] * p.unwrap(monom_i);
            coeff *= midx[monom_i] * p.m_container[monom_i];
        }
        DEBUG(midx);

        DEBUG("Coeff: " << BRY::multinom(midx) * coeff);
        p_new.m_container[p_new.wrap(sum_exp)] = BRY::multinom(midx) * coeff;
    }
    return p_new;
}