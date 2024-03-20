#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>


namespace BRY {

struct DynamicContainer {
    std::vector<bry_float_t> arr;
};

template <uint64_t SZ>
struct FixedContainer {
    std::array<bry_float_t, SZ> arr;
};

template <std::size_t DIM, class CONTAINER_T = DynamicContainer>
class Polynomial {
    public:
        Polynomial(std::size_t degree);
    public:
        template <typename ... DEGS>
        BRY_INL bry_float_t coeff(DEGS ... degrees);
    //private:
        BRY_INL std::size_t wrap(const std::array<bry_deg_t, DIM>& degrees) const;
        BRY_INL std::array<bry_deg_t, DIM> unwrap(std::size_t idx) const;
        
    private:
        std::size_t m_degree;
        CONTAINER_T m_container;
};

}

#include "impl/Polynomial_impl.hpp"