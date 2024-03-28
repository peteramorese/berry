#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>
#include <memory>

namespace BRY {

class MultiIndex {
    public:
        /// @brief Constructor
        /// @param sz Size of the index (number of individual indices)
        /// @param l1_norm Fixed sum of each index
        /// @param left Constructed at the first combination, i.e. (l1_norm, 0, ..., 0). If false, constructs at the last
        /// combination, i.e. (0, ..., 0, l1_norm)
        MultiIndex(std::size_t sz, std::size_t l1_norm, bool left = true);

        /// @brief Number of unary indices
        BRY_INL std::size_t size() const;

        /// @brief Fixed sum of each index
        BRY_INL std::size_t l1Norm() const;

        /// @brief Prefix increment. Moves the multi index along (right) by one step
        BRY_INL MultiIndex& operator++();

        /// @brief Postfix increment. Moves the multi index along (right) by one step
        BRY_INL MultiIndex operator++(int);

        /// @brief Prefix decrement. Moves the multi index backwards (left) by one step
        BRY_INL MultiIndex& operator--();

        /// @brief Postfix decrement. Moves the multi index backwards (left) by one step
        BRY_INL MultiIndex operator--(int);

        /// @brief Check if the leftmost permutation is reached
        BRY_INL bool left() const;

        /// @brief Check if the rightmost permutation is reached
        BRY_INL bool right() const;

        /// @brief Subscript operator for accessing a unary index
        /// @param d Subscript of the individual unary index
        BRY_INL std::size_t operator[](std::size_t d) const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL std::vector<std::size_t>::const_iterator begin() const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL std::vector<std::size_t>::const_iterator end() const;

    private:
        BRY_INL void updateIdx();

    private:
        std::vector<bool> m_combination;
        std::vector<std::size_t> m_idx;
        std::size_t m_l1_norm;
        bool m_left, m_right;
};

}

std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex& p);

#include "impl/MultiIndex_impl.hpp"
