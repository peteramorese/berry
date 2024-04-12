#pragma once

#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"

#include <cmath>

namespace _BRY {
    template <std::size_t DIM>
    Eigen::Tensor<BRY::bry_float_t, DIM> expandToMatchSize(const Eigen::Tensor<BRY::bry_float_t, DIM>& tensor, BRY::bry_deg_t sz) {

        #ifdef BRY_ENABLE_BOUNDS_CHECK
            ASSERT(tensor.dimension(0) < sz, "Input tensor is not smaller than desired size");
        #endif

        std::array<std::pair<BRY::bry_idx_t, BRY::bry_idx_t>, DIM> paddings;
        for (std::pair<BRY::bry_idx_t, BRY::bry_idx_t>& pads : paddings) {
            pads.first = 0;
            pads.second = sz - tensor.dimension(0);
        }

        return tensor.pad(paddings);
    }
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::Polynomial<DIM, BASIS>::Polynomial(bry_deg_t degree)
    : m_tensor(makeUniformArray<bry_deg_t, DIM>(degree + 1))
{
    m_tensor.setZero();
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::Polynomial<DIM, BASIS>::Polynomial(const Eigen::Tensor<bry_float_t, DIM>& tensor) 
    : m_tensor(tensor)
{}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::Polynomial<DIM, BASIS>::Polynomial(Eigen::Tensor<bry_float_t, DIM>&& tensor) 
    : m_tensor(std::move(tensor))
{}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_deg_t BRY::Polynomial<DIM, BASIS>::degree() const {
    return m_tensor.dimension(0) - 1;
}

template <std::size_t DIM, BRY::Basis BASIS>
template <typename ... DEGS>
BRY::bry_float_t& BRY::Polynomial<DIM, BASIS>::coeff(DEGS ... exponents) {
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_tensor(exponents...);
}

template <std::size_t DIM, BRY::Basis BASIS>
template <typename ... DEGS>
const BRY::bry_float_t& BRY::Polynomial<DIM, BASIS>::coeff(DEGS ... exponents) const {
    static_assert(is_uniform_convertible_type<bry_deg_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_deg_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_tensor(exponents...);
}

/* TODO */
//template <std::size_t DIM, BRY::Basis BASIS>
//template <typename ... FLTS>
//BRY::bry_float_t BRY::Polynomial<DIM, BASIS>::operator()(FLTS ... x) const {
//    static_assert(is_uniform_convertible_type<bry_float_t, FLTS ...>(), "All parameters passed to `operator()` must be float type (`bry_float_t`)");
//    static_assert(sizeof...(FLTS) == DIM, "Number of x parameters must match the dimension of the polynomial");
//    auto x_arr = makeArray<bry_float_t>(x ...);
//
//    //for (bry_float_t coeff)
//}


template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM, BRY::Basis::Power>& p) {
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
        if (std::abs(coeff) < BRY_OUTPUT_FMT_ZERO_THRESH) {
            iterate();
            continue;
        }

        if (!first)
            os << BRY_LOG_WHITE(" + ");
        first = false;

        os << BRY_LOG_BYELLOW(coeff);
        if constexpr (DIM != 1) {
            for (std::size_t dim = 0; dim < DIM; ++dim) {
                if (idx_arr[dim] > 0)
                    os << BRY_LOG_WHITE("(x" << dim << "^") << BRY_LOG_BGREEN(idx_arr[dim]) << BRY_LOG_WHITE(")");
            }
        } else {
            if (idx_arr[0] > 0)
                os << BRY_LOG_WHITE("x^") << BRY_LOG_BGREEN(idx_arr[0]);
        }
        iterate();
    }
    return os;
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_deg_t BRY::Polynomial<DIM, BASIS>::nMonomials() const {
    return m_tensor.size();
}

template <std::size_t DIM, BRY::Basis BASIS>
const Eigen::Tensor<BRY::bry_float_t, DIM>& BRY::Polynomial<DIM, BASIS>::tensor() const {
    return m_tensor;
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator+(BRY::bry_float_t scalar, const BRY::Polynomial<DIM, BRY::Basis::Power>& p) {
    Eigen::Tensor<BRY::bry_float_t, DIM> new_tensor = p.tensor();
    *new_tensor.data() += scalar;
    return BRY::Polynomial<DIM, BRY::Basis::Power>(std::move(new_tensor));
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator+(const BRY::Polynomial<DIM, BRY::Basis::Power>& p_1, const BRY::Polynomial<DIM, BRY::Basis::Power>& p_2) {
    const BRY::Polynomial<DIM, BRY::Basis::Power>* p_big;
    const BRY::Polynomial<DIM, BRY::Basis::Power>* p_small;
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

    Eigen::Tensor<BRY::bry_float_t, DIM> new_tensor = p_small->tensor().pad(paddings);
    new_tensor += p_big->tensor();
    BRY::Polynomial<DIM, BRY::Basis::Power> p_new(std::move(new_tensor));
    return p_new;
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator-(const BRY::Polynomial<DIM, BRY::Basis::Power>& p) {
    return BRY::Polynomial<DIM, BRY::Basis::Power>(-p.tensor());
}


template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator-(const BRY::Polynomial<DIM, BRY::Basis::Power>& p_1, const BRY::Polynomial<DIM, BRY::Basis::Power>& p_2) {
    return p_1 + -p_2;
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator*(BRY::bry_float_t scalar, const BRY::Polynomial<DIM, BRY::Basis::Power>& p) {
    return BRY::Polynomial<DIM, BRY::Basis::Power>(scalar * p.tensor());
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator*(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_float_t scalar) {
    return scalar * p;
}

/* TODO: Make this faster */
template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator*(const BRY::Polynomial<DIM, BRY::Basis::Power>& p_1, const BRY::Polynomial<DIM, BRY::Basis::Power>& p_2) {

    BRY::bry_deg_t desired_size = p_1.degree() + p_2.degree() + 2;

    Eigen::Tensor<BRY::bry_float_t, DIM> p_1_tensor_rszd = _BRY::expandToMatchSize<DIM>(p_1.tensor(), desired_size);
    Eigen::Tensor<BRY::bry_float_t, DIM> p_2_tensor_rszd = _BRY::expandToMatchSize<DIM>(p_2.tensor(), desired_size);

    std::array<BRY::bry_idx_t, DIM> dimensions;
    for (std::size_t i = 0; i < DIM; ++i)
        dimensions[i] = i;

    Eigen::Tensor<BRY::bry_complex_t, DIM> tensor_1_fft = p_1_tensor_rszd.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    Eigen::Tensor<BRY::bry_complex_t, DIM> tensor_2_fft = p_2_tensor_rszd.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    Eigen::Tensor<BRY::bry_complex_t, DIM> product_fft = tensor_1_fft * tensor_2_fft;

    Eigen::Tensor<BRY::bry_float_t, DIM> result = product_fft.template fft<Eigen::RealPart, Eigen::FFT_REVERSE>(dimensions);
    return BRY::Polynomial<DIM, BRY::Basis::Power>(std::move(result));
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator^(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_deg_t deg) {

    BRY::bry_deg_t desired_size = deg * (p.degree() + 1);

    Eigen::Tensor<BRY::bry_float_t, DIM> p_tensor_rszd = _BRY::expandToMatchSize<DIM>(p.tensor(), desired_size);

    std::array<BRY::bry_idx_t, DIM> dimensions;
    for (std::size_t i = 0; i < DIM; ++i)
        dimensions[i] = i;

    Eigen::Tensor<BRY::bry_complex_t, DIM> tensor_fft = p_tensor_rszd.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    Eigen::Tensor<BRY::bry_complex_t, DIM> exp_fft = tensor_fft.pow(static_cast<BRY::bry_float_t>(deg));

    Eigen::Tensor<BRY::bry_float_t, DIM> result = exp_fft .template fft<Eigen::RealPart, Eigen::FFT_REVERSE>(dimensions);
    return BRY::Polynomial<DIM, BRY::Basis::Power>(std::move(result));
}
