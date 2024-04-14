#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"
#include "Polynomial.h"

#include <array>

namespace BRY {

template <std::size_t DIM>
class BernsteinBasisTransform {
    public:
        /// @brief Convert a polynomial of the power basis to the Bernstein basis
        /// @param p Polynomial (power basis)
        /// @param degree_increase Elevate the degree of the Bernstein basis. 
        /// If this is 0, the degree of the Bernstein polynomial is the same as the power polynomial
        /// @return Polynomial in the Bernstein basis
        Polynomial<DIM, Basis::Bernstein> to(const Polynomial<DIM, Basis::Power>& p, bry_deg_t degree_increase = 0);

        /* TODO */
        //Polynomial<DIM, Basis::Power> from(const Polynomial<DIM, Basis::Bernstein>& p);

        

        BRY_INL bry_float_t minCoeff() const;

        /// @brief Compute the transformation matrix for power basis to Bernstein basis
        /// @param degree Degree of the power basis polynomial
        /// @param degree_increase Elevate the degree of the transformation
        /// @return Transformation matrix of elevated degree (`degree + degree_increase`)
        static Eigen::MatrixXd getTfMatrix(bry_deg_t degree, bry_deg_t degree_increase = 0);

        /// @brief Compute the inverse transformation matrix for Bernstein basis to power basis
        /// @param degree Degree of the Bernstein basis polynomial
        /// @return Transformation matrix
        static Eigen::MatrixXd getInvTfMatrix(bry_deg_t degree);
    private:
        template <typename COEFF_LAM>
        static Eigen::MatrixXd makeBigMatrix(bry_deg_t to_degree, bry_deg_t from_degree, COEFF_LAM makeCoeff);

    private:
        bry_float_t m_min_coeff;
};

}

#include "impl/BernsteinTransform_impl.hpp"