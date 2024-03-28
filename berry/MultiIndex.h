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
        /// @brief Constructor
        /// @param sz Size of the index (number of individual indices)
        /// @param l1_norm Fixed sum of each index
        /// @param begin Constructed at the first combination, i.e. (l1_norm, 0, ..., 0). If false, constructs at the last
        /// combination, i.e. (0, ..., 0, l1_norm)
        MultiIndex(std::size_t sz, std::size_t l1_norm, bool begin = true);

        BRY_INL std::size_t size() const;

        /// @brief Prefix increment. Moves the multi index along by one step
        BRY_INL MultiIndex& operator++();

        /// @brief Postfix increment. Moves the multi index along by one step
        BRY_INL MultiIndex operator++(int);

        /// @brief Prefix decrement. Moves the multi index backwards by one step
        BRY_INL MultiIndex& operator--();

        /// @brief Postfix decrement. Moves the multi index backwards by one step
        BRY_INL MultiIndex operator--(int);

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
        BRY_INL void updateIdx();

    private:
        std::vector<bool> m_combination;
        std::vector<std::size_t> m_idx;
        std::size_t m_l1_norm;
        bool m_begin, m_end;
};

}

std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex& p);

#include "impl/MultiIndex_impl.hpp"
