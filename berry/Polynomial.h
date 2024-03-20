#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>
#include <memory>


namespace BRY {

template <std::size_t DIM>
class Polynomial {
    public:
        Polynomial(bry_deg_t degree);

        template <typename ... DEGS>
        BRY_INL bry_float_t& coeff(DEGS ... exponents);

        template <typename ... DEGS>
        BRY_INL const bry_float_t& coeff(DEGS ... exponents) const;

    private:
        BRY_INL std::size_t wrap(const std::array<bry_deg_t, DIM>& exponents) const;
        BRY_INL std::array<bry_deg_t, DIM> unwrap(std::size_t idx) const;
        
    private:
        std::size_t m_degree;
        std::vector<bry_float_t> m_container;
};

}

#include "impl/Polynomial_impl.hpp"