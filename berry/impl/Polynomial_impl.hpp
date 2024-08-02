#pragma once

#include "Polynomial.h"
#include "MultiIndex.h"
#include "Operations.h"

#include "lemon/Logging.h"

#include <cmath>
#include <math.h>
#include <stdexcept>

namespace _BRY {
    template <std::size_t DIM>
    Eigen::Tensor<BRY::bry_float_t, DIM> expandToMatchSize(const Eigen::Tensor<BRY::bry_float_t, DIM>& tensor, BRY::bry_int_t sz) {
        #ifdef BRY_ENABLE_BOUNDS_CHECK
            ASSERT(tensor.dimension(0) <= sz, "Input tensor is not smaller than desired size");
        #endif

        if (tensor.dimension(0) == sz)
            return tensor;

        std::array<std::pair<BRY::bry_int_t, BRY::bry_int_t>, DIM> paddings;
        for (std::pair<BRY::bry_int_t, BRY::bry_int_t>& pads : paddings) {
            pads.first = 0;
            pads.second = sz - tensor.dimension(0);
        }

        return tensor.pad(paddings);
    }
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::Polynomial<DIM, BASIS>::Polynomial(bry_int_t degree)
    : m_tensor(makeUniformArray<bry_int_t, DIM>(degree + 1))
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
BRY::Polynomial<DIM, BASIS>::Polynomial(const Vector& vector)
{
    bry_int_t new_size = static_cast<bry_int_t>(std::pow(vector.size(), 1.0 / static_cast<bry_float_t>(DIM)));
    if (pow(new_size, DIM) < vector.size())
        new_size += 1;

    if (pow(new_size, DIM) != vector.size()) {
        ERROR("Input vector dimension mismatch");
        throw std::invalid_argument("Input vector dimension mismatch");
    }

    m_tensor = Eigen::Tensor<bry_float_t, DIM>(makeUniformArray<bry_int_t, DIM>(new_size));
    Eigen::Map<Vector> p_vec(m_tensor.data(), vector.size());
    p_vec = vector;
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_int_t BRY::Polynomial<DIM, BASIS>::degree() const {
    return m_tensor.dimension(0) - 1;
}

template <std::size_t DIM, BRY::Basis BASIS>
template <typename ... DEGS>
BRY::bry_float_t& BRY::Polynomial<DIM, BASIS>::coeff(DEGS ... exponents) {
    static_assert(is_uniform_convertible_type<bry_int_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_int_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_tensor(exponents...);
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_float_t& BRY::Polynomial<DIM, BASIS>::coeff(const std::array<bry_int_t, DIM>& exponents) {
    return m_tensor(exponents);
}

template <std::size_t DIM, BRY::Basis BASIS>
template <typename ... DEGS>
BRY::bry_float_t BRY::Polynomial<DIM, BASIS>::coeff(DEGS ... exponents) const {
    static_assert(is_uniform_convertible_type<bry_int_t, DEGS ...>(), "All parameters passed to `coeff` must be degree type (`bry_int_t`)");
    static_assert(sizeof...(DEGS) == DIM, "Number of exponents must match the dimension of the polynomial");
    return m_tensor(exponents...);
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_float_t BRY::Polynomial<DIM, BASIS>::coeff(const std::array<bry_int_t, DIM>& exponents) const {
    return m_tensor(exponents);
}

template <std::size_t DIM, BRY::Basis BASIS>
template <typename ... FLTS>
BRY::bry_float_t BRY::Polynomial<DIM, BASIS>::operator()(FLTS ... x) const {
    static_assert(is_uniform_convertible_type<bry_float_t, FLTS ...>(), "All parameters passed to `operator()` must be float type (`bry_float_t`)");
    static_assert(sizeof...(FLTS) == DIM, "Number of x parameters must match the dimension of the polynomial");
    return operator()(makeArray<bry_float_t>(x ...));
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_float_t BRY::Polynomial<DIM, BASIS>::operator()(const std::array<bry_float_t, DIM>& x) const {
    static_assert(BASIS == BRY::Basis::Power, "Evaluation of polynomials not in Power basis currently not supported");
    if constexpr (BASIS == BRY::Basis::Power) {

        // Used to store temporary sums of each x variable multiplier
        auto x_cache = makeUniformArray<bry_float_t, DIM + 1>(0.0);

        std::array<bry_int_t, DIM> deg_powers;
        for (std::size_t d = 0; d < DIM; ++d) {
            deg_powers[d] = pow(degree() + 1, d);
        }

        const bry_float_t* tensor_end_ptr = m_tensor.data() + m_tensor.size();

        for (bry_int_t i = 0; i < m_tensor.size() - 1; ++i) {
            // Set the 0'th cache spot to always be the coefficient
            x_cache[0] = *(--tensor_end_ptr);

            bry_float_t multiplier;
            bry_int_t cache_idx = DIM;
            for (; cache_idx >= 1; --cache_idx) {
                if ((i + 1) % deg_powers[cache_idx - 1] == 0) {
                    multiplier = x[cache_idx - 1];
                    break;
                }
            }

            // Set prior caches to zero
            for (bry_int_t j = 0; j < cache_idx; ++j) {
                x_cache[j + 1] += x_cache[j];
                x_cache[j] = 0.0;
            }
            
            // Multiply the current cache by the corresponding x value
            x_cache[cache_idx] *= multiplier;
        }
        x_cache[0] = *m_tensor.data();
        return std::accumulate(x_cache.begin(), x_cache.end(), 0.0);
    } else {
    }
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_float_t BRY::Polynomial<DIM, BASIS>::operator()(const Eigen::Vector<bry_float_t, DIM>& x) const {
    std::array<bry_float_t, DIM> x_arr;
    std::copy(x.begin(), x.end(), x_arr.begin());
    return operator()(x_arr);
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::Polynomial<DIM, BASIS> BRY::Polynomial<DIM, BASIS>::derivative(bry_int_t dx_idx) const {
    #ifdef BRY_ENABLE_BOUNDS_CHECK
        ASSERT(dx_idx < DIM && dx_idx >= 0, "Derivative idx out of bounds");
    #endif

    if (degree() == 0) {
        Eigen::Tensor<bry_float_t, DIM> t(m_tensor.dimensions());
        t.setZero();
        return Polynomial<DIM, BASIS>(std::move(t));
    }

    Eigen::Tensor<bry_float_t, DIM> increment_tensor = makeIncrementTensor(m_tensor.dimensions(), dx_idx, 0);

    // Multiply each coefficient in the tensor by the previous exponent (power rule)
    Eigen::Tensor<bry_float_t, DIM> power_coefficient_tensor = m_tensor * increment_tensor;

    std::array<bry_int_t, DIM> offsets = makeUniformArray<bry_int_t, DIM>(0);
    // Offset the tensor along the differential index to remove all of the 'constant' terms
    offsets[dx_idx] = 1;

    std::array<bry_int_t, DIM> extents = makeUniformArray<bry_int_t, DIM>(degree() + 1);
    // Along the differential dimension, only keep the remaining rows after deleting the first
    extents[dx_idx] -= 1;

    // Erase all the constant terms and make a tensor of dim-1 along dx_idx
    Eigen::Tensor<bry_float_t, DIM> derivative_tensor = power_coefficient_tensor.slice(offsets, extents);

    std::array<std::pair<bry_int_t, bry_int_t>, DIM> paddings;
    for (std::size_t d = 0; d < DIM; ++d) {
        if (d == dx_idx) {
            paddings[d] = std::make_pair(0, 1);
        } else {
            paddings[d] = std::make_pair(0, 0);
        }
    }

    // Add zeros on the end the tensor so that it returns to original degree (effectively)
    // shifting the coefficients over (reducing the exponent by 1)
    Eigen::Tensor<bry_float_t, DIM> derivative_tensor_orig_deg = derivative_tensor.pad(paddings);

    return Polynomial<DIM, BASIS>(std::move(derivative_tensor_orig_deg));
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::Polynomial<DIM, BASIS> BRY::Polynomial<DIM, BASIS>::liftDegree(bry_int_t raised_deg) const {
    #ifdef BRY_ENABLE_BOUNDS_CHECK
        ASSERT(raised_deg >= degree(), "Raised degree is smaller than current degree");
    #endif
    return Polynomial<DIM, BASIS>(_BRY::expandToMatchSize<DIM>(m_tensor, raised_deg + 1));
}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM, BRY::Basis::Power>& p) {
    std::array<BRY::bry_int_t, DIM> idx_arr = BRY::makeUniformArray<BRY::bry_int_t, DIM>(BRY::bry_int_t{});
    std::size_t d = 0;
    bool first = true;
    bool zero = true;

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

        zero = false;

        if (!first)
            os << LMN_LOG_WHITE(" + ");
        first = false;

        os << LMN_LOG_BYELLOW(coeff);
        if constexpr (DIM != 1) {
            for (std::size_t dim = 0; dim < DIM; ++dim) {
                if (idx_arr[dim] > 0)
                    os << LMN_LOG_WHITE("(x" << dim << "^") << LMN_LOG_BGREEN(idx_arr[dim]) << LMN_LOG_WHITE(")");
            }
        } else {
            if (idx_arr[0] > 0)
                os << LMN_LOG_WHITE("x^") << LMN_LOG_BGREEN(idx_arr[0]);
        }
        iterate();
    }

    if (zero)
        os << LMN_LOG_BYELLOW('0');
    return os;
}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM, BRY::Basis::Bernstein>& p) {
    os << "TODO";
    return os;
}

template <std::size_t DIM, BRY::Basis BASIS>
BRY::bry_int_t BRY::Polynomial<DIM, BASIS>::nMonomials() const {
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

    std::array<std::pair<BRY::bry_int_t, BRY::bry_int_t>, DIM> paddings;
    for (std::pair<BRY::bry_int_t, BRY::bry_int_t>& pads : paddings) {
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

    BRY::bry_int_t desired_size = p_1.degree() + p_2.degree() + 1;

    Eigen::Tensor<BRY::bry_float_t, DIM> p_1_tensor_rszd = _BRY::expandToMatchSize<DIM>(p_1.tensor(), desired_size);
    Eigen::Tensor<BRY::bry_float_t, DIM> p_2_tensor_rszd = _BRY::expandToMatchSize<DIM>(p_2.tensor(), desired_size);

    std::array<BRY::bry_int_t, DIM> dimensions;
    for (std::size_t i = 0; i < DIM; ++i)
        dimensions[i] = i;

    Eigen::Tensor<BRY::bry_complex_t, DIM> tensor_1_fft = p_1_tensor_rszd.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    Eigen::Tensor<BRY::bry_complex_t, DIM> tensor_2_fft = p_2_tensor_rszd.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    Eigen::Tensor<BRY::bry_complex_t, DIM> product_fft = tensor_1_fft * tensor_2_fft;

    Eigen::Tensor<BRY::bry_float_t, DIM> result = product_fft.template fft<Eigen::RealPart, Eigen::FFT_REVERSE>(dimensions);
    return BRY::Polynomial<DIM, BRY::Basis::Power>(std::move(result));
}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator^(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_int_t exp) {

    if (exp == 0) {
        Eigen::Tensor<BRY::bry_float_t, DIM> scalar_t(makeUniformArray<BRY::bry_int_t, DIM>(1));
        *scalar_t.data() = 1;
        return BRY::Polynomial<DIM, BRY::Basis::Power>(scalar_t);
    }

    BRY::bry_int_t desired_size = exp * p.degree() + 1;

    Eigen::Tensor<BRY::bry_float_t, DIM> p_tensor_rszd = _BRY::expandToMatchSize<DIM>(p.tensor(), desired_size);

    std::array<BRY::bry_int_t, DIM> dimensions;
    for (std::size_t i = 0; i < DIM; ++i)
        dimensions[i] = i;

    Eigen::Tensor<BRY::bry_complex_t, DIM> tensor_fft = p_tensor_rszd.template fft<Eigen::BothParts, Eigen::FFT_FORWARD>(dimensions);
    Eigen::Tensor<BRY::bry_complex_t, DIM> exp_fft = tensor_fft.pow(static_cast<BRY::bry_float_t>(exp));

    Eigen::Tensor<BRY::bry_float_t, DIM> result = exp_fft .template fft<Eigen::RealPart, Eigen::FFT_REVERSE>(dimensions);
    return BRY::Polynomial<DIM, BRY::Basis::Power>(std::move(result));
}

template <std::size_t DIM, BRY::Basis FROM_BASIS, BRY::Basis TO_BASIS>
BRY::Polynomial<DIM, TO_BASIS> BRY::transform(const Polynomial<DIM, FROM_BASIS>& p, const Matrix& transform_matrix) {
    bry_int_t new_size = static_cast<bry_int_t>(std::pow(transform_matrix.rows(), 1.0 / static_cast<bry_float_t>(DIM)));
    if (pow(new_size, DIM) < transform_matrix.rows())
        new_size += 1;

    Eigen::Tensor<bry_float_t, DIM> tensor(makeUniformArray<bry_int_t, DIM>(new_size));
    tensor.setZero();

    Eigen::Map<const Vector> p_vec(p.tensor().data(), p.nMonomials());
    Eigen::Map<Vector> p_vec_tf(tensor.data(), transform_matrix.rows());

    p_vec_tf = transform_matrix * p_vec;

    return BRY::Polynomial<DIM, TO_BASIS>(std::move(tensor));
}
