#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>
#include <memory>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

/* Threshold for which small coefficients will be ignored when printing
    NOTE: Does not affect numerical calculations
*/
#define BRY_OUTPUT_FMT_ZERO_THRESH 0.00000001

/// Forward declarations
namespace BRY {

template <std::size_t DIM>
class Polynomial;

}

template <std::size_t DIM>
BRY::Polynomial<DIM> operator+(BRY::bry_float_t scalar, const BRY::Polynomial<DIM>& p);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator+(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator-(const BRY::Polynomial<DIM>& p);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator-(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(BRY::bry_float_t scalar, const BRY::Polynomial<DIM>& p);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(const BRY::Polynomial<DIM>& p, BRY::bry_float_t scalar);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator*(const BRY::Polynomial<DIM>& p_1, const BRY::Polynomial<DIM>& p_2);

template <std::size_t DIM>
BRY::Polynomial<DIM> operator^(const BRY::Polynomial<DIM>& p, BRY::bry_deg_t deg);

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM>& p);

namespace BRY {

template <std::size_t DIM>
class Polynomial {
    public:
        Polynomial(bry_deg_t degree);

        BRY_INL std::size_t degree() const;

        template <typename ... DEGS>
        BRY_INL bry_float_t& coeff(DEGS ... exponents);

        template <typename ... DEGS>
        BRY_INL const bry_float_t& coeff(DEGS ... exponents) const;

        // Operators

        /* TODO */
        template <typename ... FLTS>
        bry_float_t operator()(FLTS ... x) const;

        friend std::ostream& operator<<<DIM>(std::ostream& os, const Polynomial& p);

    private:
        Polynomial(const Eigen::Tensor<bry_float_t, DIM>& tensor);
        Polynomial(Eigen::Tensor<bry_float_t, DIM>&& tensor);

    private:
        Eigen::Tensor<bry_float_t, DIM> m_tensor;

    private:
        //template <std::size_t _DIM>
        friend Polynomial operator+<DIM>(bry_float_t scalar, const Polynomial& p);
        friend Polynomial operator+<DIM>(const Polynomial& p_1, const Polynomial& p_2);
        friend Polynomial operator-<DIM>(const Polynomial& p);
        friend Polynomial operator-<DIM>(const Polynomial& p_1, const Polynomial& p_2);
        friend Polynomial operator*<DIM>(bry_float_t scalar, const BRY::Polynomial<DIM>& p);
        friend Polynomial operator*<DIM>(const Polynomial& p_1, const Polynomial& p_2);
        friend Polynomial operator^<DIM>(const Polynomial& p, bry_deg_t deg);
};

}

#include "impl/Polynomial_impl.hpp"