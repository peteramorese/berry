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

namespace BRY {

/// @brief Different supported polynomial bases
enum class Basis {
    Power,
    Bernstein
};

}

/* Forward Declarations */
namespace BRY {

template <std::size_t DIM, BRY::Basis BASIS>
class Polynomial;

}

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator+(BRY::bry_float_t scalar, const BRY::Polynomial<DIM, BRY::Basis::Power>& p);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator+(const BRY::Polynomial<DIM, BRY::Basis::Power>& p_1, const BRY::Polynomial<DIM, BRY::Basis::Power>& p_2);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator-(const BRY::Polynomial<DIM, BRY::Basis::Power>& p);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator-(const BRY::Polynomial<DIM, BRY::Basis::Power>& p_1, const BRY::Polynomial<DIM, BRY::Basis::Power>& p_2);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator*(BRY::bry_float_t scalar, const BRY::Polynomial<DIM, BRY::Basis::Power>& p);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator*(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_float_t scalar);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator*(const BRY::Polynomial<DIM, BRY::Basis::Power>& p_1, const BRY::Polynomial<DIM, BRY::Basis::Power>& p_2);

template <std::size_t DIM>
BRY::Polynomial<DIM, BRY::Basis::Power> operator^(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_deg_t deg);

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM, BRY::Basis::Power>& p);

namespace BRY {

template <std::size_t DIM, Basis BASIS = Basis::Power>
class Polynomial {
    public:
        /// @brief Construct polynomial of known (max) degree
        /// @param degree Maximum exponent of a given variable
        Polynomial(bry_deg_t degree);

        /// @brief Construct polynomial given a tensor
        /// @param tensor Tensor of coefficients 
        Polynomial(const Eigen::Tensor<bry_float_t, DIM>& tensor);

        /// @brief Construct polynomial moving a tensor
        /// @param tensor Tensor of coefficients
        Polynomial(Eigen::Tensor<bry_float_t, DIM>&& tensor);

        /// @brief Get the degree
        /// @return Degree
        BRY_INL bry_deg_t degree() const;

        /// @brief Access a specific coefficient of a term. Usage (power basis): coeff(1, 0, 3) returns the coefficient of the term (x0)(x2^3)
        /// @tparam ...DEGS 
        /// @param ...exponents Exponents in order of variables
        /// @return Reference to mutable value
        template <typename ... DEGS>
        BRY_INL bry_float_t& coeff(DEGS ... exponents);

        /// @brief Access a specific coefficient of a term. Usage (power basis): coeff(1, 0, 3) returns the coefficient of the term (x0)(x2^3)
        /// @tparam ...DEGS 
        /// @param ...exponents Exponents in order of variables
        /// @return Reference to imutable value
        template <typename ... DEGS>
        BRY_INL const bry_float_t& coeff(DEGS ... exponents) const;

        // Operators

        /* TODO */
        template <typename ... FLTS>
        bry_float_t operator()(FLTS ... x) const;

        friend std::ostream& operator<<<DIM>(std::ostream& os, const Polynomial& p);

        BRY_INL const Eigen::Tensor<bry_float_t, DIM>& tensor() const;

    private:

    private:
        Eigen::Tensor<bry_float_t, DIM> m_tensor;

};

}

#include "impl/Polynomial_impl.hpp"