#pragma once

#include "Options.h"
#include "Types.h"
#include "Polynomial.h"

#include <array>

namespace BRY {

template <std::size_t DIM>
class BernsteinBasisTransform {
    public:
        /// @brief Compute the transformation matrix for power basis to Bernstein basis
        /// @param degree Degree of the power basis polynomial
        /// @param degree_increase Elevate the degree of the transformation
        /// @return Transformation matrix of elevated degree (`degree + degree_increase`)
        static Matrix pwrToBernMatrix(bry_int_t degree, bry_int_t degree_increase = 0);

        /// @brief Compute the inverse transformation matrix for Bernstein basis to power basis
        /// @param degree Degree of the Bernstein basis polynomial
        /// @return Transformation matrix
        static Matrix bernToPwrMatrix(bry_int_t degree);

        /// @brief Compute the lower bound of a polynomial on the unit interval in the Bernstein basis
        /// @param p Polynomial in the Bernstein basis
        /// @return Lower bound (smallest coefficient), flag if the vertex condition is met (true lower bound achieved)
        static std::pair<bry_float_t, bool> infBound(const BRY::Polynomial<DIM, BRY::Basis::Bernstein>& p);

        /// @brief Compute the lower bound of a polynomial on the unit interval in the Bernstein basis, track the idx of the mim coeff
        /// @param p Polynomial in the Bernstein basis
        /// @param coefficient_idx Edit in place (return) the idx of the min coefficient
        /// @return Lower bound (smallest coefficient), flag if the vertex condition is met (true lower bound achieved)
        static std::pair<bry_float_t, bool> infBound(const BRY::Polynomial<DIM, BRY::Basis::Bernstein>& p, std::array<bry_int_t, DIM>& coefficient_idx);

        /// @brief Compute the difference between an upper and lower bound on the infemum of a polynomial
        /// @param p Polynomial in the power basis
        /// @param degree_increase Elevated degree of Bernstein transformation
        /// @return Gap between inf lower and upper bound
        static bry_float_t infBoundGap(const BRY::Polynomial<DIM, BRY::Basis::Power>& p, bool vertex_condition = false, bry_int_t degree_increase = 0);

        static Eigen::Vector<bry_float_t, DIM> ctrlPtOnUnitBox(const std::array<bry_int_t, DIM>& coefficient_idx, bry_int_t bernstein_p_deg);

    private:
        template <typename COEFF_LAM>
        static Matrix makeBigMatrix(bry_int_t to_degree, bry_int_t from_degree, COEFF_LAM makeCoeff);

    private:
        bry_float_t m_min_coeff;
};

}

#include "impl/BernsteinTransform_impl.hpp"