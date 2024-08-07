#pragma once

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
BRY::Polynomial<DIM, BRY::Basis::Power> operator^(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_int_t exp);

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM, BRY::Basis::Power>& p);

template <std::size_t DIM>
std::ostream& operator<<(std::ostream& os, const BRY::Polynomial<DIM, BRY::Basis::Bernstein>& p);

namespace BRY {

template <std::size_t DIM, Basis BASIS = Basis::Power>
class Polynomial {
    public:
        /// @brief Construct polynomial of known (max) degree
        /// @param degree Maximum exponent of a given variable
        Polynomial(bry_int_t degree);

        /// @brief Construct polynomial given a tensor
        /// @param tensor Tensor of coefficients 
        Polynomial(const Eigen::Tensor<bry_float_t, DIM>& tensor);

        /// @brief Construct polynomial moving a tensor
        /// @param tensor Tensor of coefficients
        Polynomial(Eigen::Tensor<bry_float_t, DIM>&& tensor);

        /// @brief Construct polynomial given a multiindex-organized vector of coefficients
        /// @param tensor Tensor of coefficients 
        Polynomial(const Vector& vector);

        /// @brief Get the degree
        /// @return Degree
        BRY_INL bry_int_t degree() const;

        /// @brief Access a specific coefficient of a term. Usage (power basis): coeff(1, 0, 3) returns the coefficient of the term (x0)(x2^3)
        /// @tparam ...DEGS 
        /// @param ...exponents Exponents in order of variables
        /// @return Reference to mutable value
        template <typename ... DEGS>
        BRY_INL bry_float_t& coeff(DEGS ... exponents);
        BRY_INL bry_float_t& coeff(const std::array<bry_int_t, DIM>& exponents);

        /// @brief Access a specific coefficient of a term. Usage (power basis): coeff(1, 0, 3) returns the coefficient of the term (x0)(x2^3)
        /// @tparam ...DEGS 
        /// @param ...exponents Exponents in order of variables
        /// @return Reference to imutable value
        template <typename ... DEGS>
        BRY_INL bry_float_t coeff(DEGS ... exponents) const;
        BRY_INL bry_float_t coeff(const std::array<bry_int_t, DIM>& exponents) const;

        /// @brief Evaluate the polynomial for given x vector
        /// @tparam ...FLTS 
        /// @param ...x `x` values
        /// @return Scalar 
        template <typename ... FLTS>
        BRY_INL bry_float_t operator()(FLTS ... x) const;
        bry_float_t operator()(const std::array<bry_float_t, DIM>& x) const;
        BRY_INL bry_float_t operator()(const Eigen::Vector<bry_float_t, DIM>& x) const;

        /// @brief Compute the (partial) derivative of the polynomial with respect to a given dimension
        /// @param dx_idx Dimension to take the partial derivative with respect to
        /// @return Derivative polynomial (with the same degree)
        Polynomial<DIM, BASIS> derivative(bry_int_t dx_idx) const;

        /// @brief Create a copy with a raised degree by padding the higher order terms as zero-coefficients
        /// @param raised_deg New degree (must be larger than `degree()`)
        /// @return Raised degree polynomial
        Polynomial<DIM, BASIS> liftDegree(bry_int_t raised_deg) const;

        /// @brief Get the Number of monomials
        bry_int_t nMonomials() const;

        /// @brief Access the underlying tensor
        /// @return Read-only tensor access
        BRY_INL const Eigen::Tensor<bry_float_t, DIM>& tensor() const;

        friend std::ostream& operator<<<DIM>(std::ostream& os, const Polynomial& p);

    private:
        Eigen::Tensor<bry_float_t, DIM> m_tensor;

};

/// @brief Linearly transform the coefficients of a polynomial using a transformation matrix
/// @tparam FROM_BASIS Basis of existing polynomial
/// @tparam TO_BASIS Basis of returned polynomial
/// @param p Polynomial
/// @param transform_matrix Transformation in vectorized form
/// @return Transformed polynomial
template <std::size_t DIM, Basis FROM_BASIS, Basis TO_BASIS = Basis::Power>
BRY::Polynomial<DIM, TO_BASIS> transform(const BRY::Polynomial<DIM, FROM_BASIS>& p, const Matrix& transform_matrix);

}

#include "impl/Polynomial_impl.hpp"