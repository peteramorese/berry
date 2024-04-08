#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>
#include <memory>

namespace BRY {

/// @brief Increment the multi index through all indices that have a fixed L-1 norm
class FixedNormIncrementer {
    public:
        FixedNormIncrementer(std::size_t sz, bry_idx_t index_constraint, bool left, std::vector<bry_idx_t>& initial_idx);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);

    private:
        BRY_INL void updateIdx(std::vector<bry_idx_t>& current_idx) const;

    private:
        bry_idx_t m_l1_norm;
        std::vector<bool> m_combination;
};

/// @brief Increment the multi index through all indices that have a L-infinity norm upper bound
class ExhaustiveIncrementer {
    public:
        ExhaustiveIncrementer(std::size_t sz, bry_idx_t index_constraint, bool left, std::vector<bry_idx_t>& initial_idx);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);
    private:
        bry_idx_t m_linfty_norm;
};

/// @brief Increment the multi index through all indices that have a L-infinity norm upper bound
class BoundedExhaustiveIncrementer {
    public:
        BoundedExhaustiveIncrementer(std::size_t sz, const std::vector<bry_idx_t>& index_bounds, bool left);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);
    private:
        bry_idx_t m_linfty_norm;
};

template <class INCREMENTER = ExhaustiveIncrementer>
class MultiIndex {
    public:
        /// @brief Automatic incrementer constructor (Supported for FixedNormIncrementer and ExhaustiveIncrementer)
        /// @param sz Size of the index (number of individual indices)
        /// @param index_constraint Norm constraint passed to the incrementer
        /// @param left Constructed at the first combination if true, last combination otherwise
        /// combination, i.e. (0, ..., 0, l1_norm)
        MultiIndex(std::size_t sz, bry_idx_t index_constraint, bool left = true);

        /// @brief Generic incrementer constructor (Supported for any custom incrementer)
        /// @param incrementer Incrementer object 
        /// @param initial_idx Initial index
        MultiIndex(INCREMENTER&& incrementer, std::vector<bry_idx_t>&& initial_idx);

        /// @brief Number of unary indices
        BRY_INL std::size_t size() const;

        /// @brief Max norm constraint passed to incrementer
        BRY_INL INCREMENTER& incrementer();
        BRY_INL const INCREMENTER& incrementer() const;

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
        BRY_INL bry_idx_t operator[](std::size_t d) const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL std::vector<bry_idx_t>::const_iterator begin() const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL std::vector<bry_idx_t>::const_iterator end() const;

    private:
        std::vector<bry_idx_t> m_idx;
        INCREMENTER m_incrementer;
        bool m_left, m_right;
};

}

template <class INCREMENTER>
std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex<INCREMENTER>& p);

#include "impl/MultiIndex_impl.hpp"
