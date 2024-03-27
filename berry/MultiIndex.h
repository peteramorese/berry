#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>
#include <memory>
#include <iostream>

namespace BRY {

class MultiIndex {
    public:
        MultiIndex(std::size_t sz, std::size_t l1_norm);

        /// @brief Prefix increment. Moves the multi index along by one step
        BRY_INL MultiIndex& operator++();

        /// @brief Postfix increment. Moves the multi index along by one step
        BRY_INL MultiIndex operator++(int);

        /// @brief Check if the index is done incrementing
        /// @return `true` if done, `false` otherwise
        BRY_INL bool end() const;

        /// @brief Subscript operator for accessing a unary index
        /// @param i Subscript of the individual unary index
        BRY_INL std::size_t operator[](std::size_t i) const;
    private:
        std::vector<std::size_t> m_midx;
        std::size_t m_l1_norm;
        std::size_t m_i;
};

}

#include "impl/MultiIndex_impl.hpp"
