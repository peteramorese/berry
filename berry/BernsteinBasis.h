#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"
#include "Polynomial.h"

#include <array>

namespace BRY {

template <std::size_t DIM>
class BernsteinBasis {
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

    private:
        Eigen::Tensor<bry_float_t, DIM> getDenominatorTensor()

    private:
        bry_float_t m_min_coeff;
};

}

#include "impl/Bernstein_impl.hpp"