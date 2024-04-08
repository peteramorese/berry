#pragma once

#include "Logging.h"
#include "Options.h"
#include "Types.h"

#include <vector>
#include <array>
#include <memory>

namespace BRY {

/// @brief Increment the multi index through all indices that have a L-infinity norm upper bound
class ExhaustiveIncrementer {
    public:
        ExhaustiveIncrementer(std::vector<bry_idx_t>& initial_idx, bool first,
            std::size_t sz, bry_idx_t index_constraint);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);
    private:
        bry_idx_t m_linfty_norm;
};

/// @brief Increment the multi index through all indices that have a L-infinity norm upper bound (keeps track of the wrapped index)
class ExhaustiveIncrementerWrap {
    public:
        ExhaustiveIncrementerWrap(std::vector<bry_idx_t>& initial_idx, bool first,
            std::size_t sz, bry_idx_t index_constraint);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);

        /// @brief This incrementer keeps track of the wrapped index (equivalent flattened 1D index) given the index constraint
        /// @return Wrapped index
        BRY_INL bry_idx_t wrappedIdx() const;
    private:
        bry_idx_t m_linfty_norm;
        bry_idx_t m_wrapped_idx;
};

/// @brief Increment the multi index through all indices that have a fixed L-1 norm
class FixedNormIncrementer {
    public:
        FixedNormIncrementer(std::vector<bry_idx_t>& initial_idx, bool first,
            std::size_t sz, bry_idx_t index_constraint);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);

    private:
        BRY_INL void updateIdx(std::vector<bry_idx_t>& current_idx) const;

    private:
        bry_idx_t m_l1_norm;
        std::vector<bool> m_combination;
};

/// @brief Increment the multi index through all indices that are less than a multi-index bound (keeps track of the wrapped index)
class BoundedExhaustiveIncrementerWrap {
    public:
        BoundedExhaustiveIncrementerWrap(std::vector<bry_idx_t>& initial_idx, bool first, 
            const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint);
        BRY_INL const std::vector<bry_idx_t>& indexConstraint() const;
        BRY_INL bool increment(std::vector<bry_idx_t>& current_idx);
        BRY_INL bool decrement(std::vector<bry_idx_t>& current_idx);

        /// @brief This incrementer keeps track of the wrapped index (equivalent flattened 1D index) given the index constraint
        /// @return Wrapped index
        BRY_INL bry_idx_t wrappedIdx() const;
    private:
        std::vector<bry_idx_t> m_index_bounds;
        bry_idx_t m_linfty_norm;
        bry_idx_t m_wrapped_idx;
};


template <class INCREMENTER = ExhaustiveIncrementer>
class MultiIndex {
    public:
        template <typename ... ARGS>
        MultiIndex(bool first, ARGS&&... args);

        /// @brief Generic incrementer constructor (Supported for any custom incrementer)
        /// @param incrementer Incrementer object 
        /// @param initial_idx Initial index
        //MultiIndex(INCREMENTER&& incrementer, std::vector<bry_idx_t>&& initial_idx, bool left = true);

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

        /// @brief Check if the first permutation is reached
        BRY_INL bool first() const;

        /// @brief Check if the last permutation is reached
        BRY_INL bool last() const;

        /// @brief Subscript operator for accessing a unary index
        /// @param d Subscript of the individual unary index
        BRY_INL bry_idx_t operator[](std::size_t d) const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL std::vector<bry_idx_t>::const_iterator begin() const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL std::vector<bry_idx_t>::const_iterator end() const;

        /// @brief Access the multi index as a vector
        BRY_INL const std::vector<bry_idx_t>& getIdxVector() const;
    private:
        std::vector<bry_idx_t> m_idx;
        INCREMENTER m_incrementer;
        bool m_first, m_last;
};

/* Convenience methods for creating MultiIndices easily (put in the for loop, don't pollute the stack!) */
BRY_INL static MultiIndex<ExhaustiveIncrementer> mIdx(std::size_t sz, bry_idx_t index_constraint);
BRY_INL static MultiIndex<ExhaustiveIncrementer> rmIdx(std::size_t sz, bry_idx_t index_constraint);
BRY_INL static MultiIndex<ExhaustiveIncrementerWrap> mIdxW(std::size_t sz, bry_idx_t index_constraint);
BRY_INL static MultiIndex<ExhaustiveIncrementerWrap> rmIdxW(std::size_t sz, bry_idx_t index_constraint);
BRY_INL static MultiIndex<FixedNormIncrementer> mIdxFN(std::size_t sz, bry_idx_t index_constraint);
BRY_INL static MultiIndex<FixedNormIncrementer> rmIdxFN(std::size_t sz, bry_idx_t index_constraint);
BRY_INL static MultiIndex<BoundedExhaustiveIncrementerWrap> mIdxBEW(const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint);
BRY_INL static MultiIndex<BoundedExhaustiveIncrementerWrap> rmIdxBEW(const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint);

}

template <class INCREMENTER>
std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex<INCREMENTER>& p);

#include "impl/MultiIndex_impl.hpp"
