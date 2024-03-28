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

        BRY_INL std::size_t size() const;

        /// @brief Prefix increment. Moves the multi index along by one step
        BRY_INL MultiIndex& operator++();

        /// @brief Postfix increment. Moves the multi index along by one step
        BRY_INL MultiIndex operator++(int);

        /// @brief Check if the initial permutation is reached
        /// @return `true` if done, `false` otherwise
        BRY_INL bool begin() const;

        /// @brief Check if the final permutation is reached
        /// @return `true` if done, `false` otherwise
        BRY_INL bool end() const;

        /// @brief Subscript operator for accessing a unary index
        /// @param d Subscript of the individual unary index
        BRY_INL std::size_t operator[](std::size_t d) const;

    private:
        std::vector<bool> m_midx;
        std::size_t m_l1_norm;
        bool m_begin, m_end;
};

}

std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex& p);

#include "impl/MultiIndex_impl.hpp"
