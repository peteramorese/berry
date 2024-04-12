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
        ExhaustiveIncrementer(bry_idx_t*& _initial_idx, std::size_t sz, bool first,
            bry_idx_t index_constraint);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(bry_idx_t* current_idx, std::size_t sz);
        BRY_INL bool decrement(bry_idx_t* current_idx, std::size_t sz);
    private:
        bry_idx_t m_linfty_norm;
};

/// @brief Increment the multi index through all indices that have a L-infinity norm upper bound (keeps track of the wrapped index)
class ExhaustiveIncrementerWrap {
    public:
        ExhaustiveIncrementerWrap(bry_idx_t*& _initial_idx, std::size_t sz, bool first,
            bry_idx_t index_constraint);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(bry_idx_t* current_idx, std::size_t sz);
        BRY_INL bool decrement(bry_idx_t* current_idx, std::size_t sz);

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
        FixedNormIncrementer(bry_idx_t*& _initial_idx, std::size_t sz, bool first,
            bry_idx_t index_constraint);
        BRY_INL bry_idx_t indexConstraint() const;
        BRY_INL bool increment(bry_idx_t* current_idx, std::size_t sz);
        BRY_INL bool decrement(bry_idx_t* current_idx, std::size_t sz);

    private:
        BRY_INL void updateIdx(bry_idx_t* current_idx, std::size_t sz) const;

    private:
        bry_idx_t m_l1_norm;
        std::vector<bool> m_combination;
};

/// @brief Increment the multi index through all indices that are less than a multi-index bound (keeps track of the wrapped index)
class BoundedExhaustiveIncrementerWrap {
    public:
        BoundedExhaustiveIncrementerWrap(bry_idx_t*& _initial_idx, std::size_t sz, bool first, 
            const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint);
        BRY_INL const std::vector<bry_idx_t>& indexConstraint() const;
        BRY_INL bool increment(bry_idx_t* current_idx, std::size_t sz);
        BRY_INL bool decrement(bry_idx_t* current_idx, std::size_t sz);

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
        /// @brief Constructor with internal array storage
        /// @param first If true, start at the first index, otherwise starts at the last index
        /// @param ...args Parameters to constructor of specific incrementer
        template <typename ... ARGS>
        MultiIndex(std::size_t sz, bool first, ARGS&&... args);

        /// @brief Constructor with external array storage
        /// @param array Pointer to external contiguous array where multi index will be edited in place
        /// @param sz Size of external array
        /// @param first If true, start at the first index, otherwise starts at the last index
        /// @param ...args Parameters to constructor of specific incrementer
        template <typename ... ARGS>
        MultiIndex(bry_idx_t* array, std::size_t sz, bool first, ARGS&&... args);

        ~MultiIndex();

        /// @brief Generic incrementer constructor (Supported for any custom incrementer)
        /// @param incrementer Incrementer object 
        /// @param initial_idx Initial index
        //MultiIndex(INCREMENTER&& incrementer, std::vector<bry_idx_t>&& initial_idx, bool left = true);

        /// @brief Number of unary indices
        BRY_INL std::size_t size() const;

        /// @brief Max norm constraint passed to incrementer
        BRY_INL INCREMENTER& inc();
        BRY_INL const INCREMENTER& inc() const;

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
        BRY_INL const bry_idx_t* begin() const;

        /// @brief Vector iterator access of the unary indices
        BRY_INL const bry_idx_t* end() const;

    private:
        bry_idx_t* m_idx;
        std::size_t m_sz;
        INCREMENTER m_incrementer;
        bool m_first, m_last;
        bool m_external_arr;
};

/* Convenience methods for creating MultiIndices easily */
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
