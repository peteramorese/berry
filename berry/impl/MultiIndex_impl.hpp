#pragma once

#include "MultiIndex.h"

#include <cmath>
#include <algorithm>

/* Exhaustive Norm Incrementer */

BRY::ExhaustiveIncrementer::ExhaustiveIncrementer(std::vector<bry_idx_t>& initial_idx, bool first, std::size_t sz, bry_idx_t index_constraint)
    : m_linfty_norm(index_constraint)
{
    if (!first) {
        for (bry_idx_t& i : initial_idx)
            i = index_constraint - 1;
    } else {
        initial_idx.resize(sz, 0);
    }
}

BRY::bry_idx_t BRY::ExhaustiveIncrementer::indexConstraint() const {
    return m_linfty_norm;
}


bool BRY::ExhaustiveIncrementer::increment(std::vector<bry_idx_t>& current_idx) {
    ++current_idx[0];
    for (std::size_t i = 0; i < current_idx.size() - 1; ++i) {
        if (current_idx[i] >= m_linfty_norm) {
            current_idx[i] = 0;
            ++current_idx[i + 1];
        }
    }
    return current_idx.back() < m_linfty_norm;
}

bool BRY::ExhaustiveIncrementer::decrement(std::vector<bry_idx_t>& current_idx) {
    --current_idx[0];
    for (std::size_t i = 0; i < current_idx.size() - 1; ++i) {
        if (current_idx[i] < 0) {
            current_idx[i] = m_linfty_norm - 1;
            --current_idx[i + 1];
        }
    }
    return current_idx.back() >= 0;
}

/* Exhaustive Norm Incrementer */

BRY::ExhaustiveIncrementerWrap::ExhaustiveIncrementerWrap(std::vector<bry_idx_t>& initial_idx, bool first, std::size_t sz, bry_idx_t index_constraint)
    : m_linfty_norm(index_constraint)
{
    if (!first) {
        for (bry_idx_t& i : initial_idx)
            i = index_constraint - 1;
        m_wrapped_idx = std::pow(m_linfty_norm, sz) - 1;
    } else {
        initial_idx.resize(sz, 0);
        m_wrapped_idx = 0;
    }
}

BRY::bry_idx_t BRY::ExhaustiveIncrementerWrap::indexConstraint() const {
    return m_linfty_norm;
}


bool BRY::ExhaustiveIncrementerWrap::increment(std::vector<bry_idx_t>& current_idx) {
    ++current_idx[0];
    ++m_wrapped_idx;
    for (std::size_t i = 0; i < current_idx.size() - 1; ++i) {
        if (current_idx[i] >= m_linfty_norm) {
            current_idx[i] = 0;
            ++current_idx[i + 1];
        }
    }
    return current_idx.back() < m_linfty_norm;
}

bool BRY::ExhaustiveIncrementerWrap::decrement(std::vector<bry_idx_t>& current_idx) {
    --current_idx[0];
    --m_wrapped_idx;
    for (std::size_t i = 0; i < current_idx.size() - 1; ++i) {
        if (current_idx[i] < 0) {
            current_idx[i] = m_linfty_norm - 1;
            --current_idx[i + 1];
        }
    }
    return current_idx.back() >= 0;
}

BRY::bry_idx_t BRY::ExhaustiveIncrementerWrap::wrappedIdx() const {
    return m_wrapped_idx;
}

/* Fixed Norm Incrementer */

BRY::FixedNormIncrementer::FixedNormIncrementer(std::vector<bry_idx_t>& initial_idx, bool first, std::size_t sz, bry_idx_t index_constraint)
    : m_combination(index_constraint + sz - 2, false)
    , m_l1_norm(index_constraint)
{
    if (first) {
        for (std::size_t i = 0; i < index_constraint - 1; ++i) {
            m_combination[i] = true;
        }
        initial_idx.front() = index_constraint - 1;
    } else {
        auto it = m_combination.end();
        for (std::size_t i = 0; i < index_constraint - 1; ++i) {
            *(--it) = true;
        }
        initial_idx.back() = index_constraint - 1;
    }
}

BRY::bry_idx_t BRY::FixedNormIncrementer::indexConstraint() const {
    return m_l1_norm;
}

bool BRY::FixedNormIncrementer::increment(std::vector<bry_idx_t>& current_idx) {
    if (std::next_permutation(m_combination.begin(), m_combination.end(), std::greater())) {
        updateIdx(current_idx);
        return true;
    } else {
        std::prev_permutation(m_combination.begin(), m_combination.end(), std::greater());
        return false;
    }
}

bool BRY::FixedNormIncrementer::decrement(std::vector<bry_idx_t>& current_idx) {
    if (std::prev_permutation(m_combination.begin(), m_combination.end(), std::greater())) {
        updateIdx(current_idx);
        return true;
    } else {
        std::next_permutation(m_combination.begin(), m_combination.end(), std::greater());
        return false;
    }
}

void BRY::FixedNormIncrementer::updateIdx(std::vector<bry_idx_t>& current_idx) const {
    std::size_t d_i = 0;
    std::size_t sum = 0;
    for (bool bit : m_combination) {
        if (bit) {
            ++sum;
        } else {
            current_idx[d_i++] = sum;
            sum = 0;
        }
    }
    current_idx.back() = sum;
}

/* Bounded Exhaustive Incrementer Wrap */

BRY::BoundedExhaustiveIncrementerWrap::BoundedExhaustiveIncrementerWrap(std::vector<bry_idx_t>& initial_idx, bool first, const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint) 
    : m_index_bounds(index_bounds)
    , m_linfty_norm(index_constraint)
{
    if (first) {
        m_wrapped_idx = 0;
        initial_idx.resize(index_bounds.size(), 0);
    } else {
        initial_idx.resize(index_bounds.size());
        m_wrapped_idx = 0;
        for (std::size_t i = 0; i < index_bounds.size(); ++i) {
            initial_idx[i] = index_bounds[i] - 1;
            m_wrapped_idx += std::pow(m_linfty_norm, i) * (index_bounds[i] - 1);
        }
    }
}

const std::vector<BRY::bry_idx_t>& BRY::BoundedExhaustiveIncrementerWrap::indexConstraint() const {
    return m_index_bounds;
}

bool BRY::BoundedExhaustiveIncrementerWrap::increment(std::vector<bry_idx_t>& current_idx) {
    ++current_idx[0];
    ++m_wrapped_idx;
    for (std::size_t i = 0; i < current_idx.size() - 1; ++i) {
        if (current_idx[i] >= m_index_bounds[i]) {
            current_idx[i] = 0;
            ++current_idx[i + 1];
            m_wrapped_idx += std::pow(m_linfty_norm, i) * (m_linfty_norm - m_index_bounds[i]);
        }
    }
    return current_idx.back() < m_index_bounds.back();
}

bool BRY::BoundedExhaustiveIncrementerWrap::decrement(std::vector<bry_idx_t>& current_idx) {
    --current_idx[0];
    --m_wrapped_idx;
    for (std::size_t i = 0; i < current_idx.size() - 1; ++i) {
        if (current_idx[i] < 0) {
            current_idx[i] = m_index_bounds[i] - 1;
            --current_idx[i + 1];
            m_wrapped_idx -= std::pow(m_linfty_norm, i) * (m_linfty_norm - m_index_bounds[i]);
        }
    }
    return current_idx.back() >= 0;
}

BRY::bry_idx_t BRY::BoundedExhaustiveIncrementerWrap::wrappedIdx() const {
    return m_wrapped_idx;
}

/* MultiIndex */

template <class INCREMENTER>
template <typename ... ARGS>
BRY::MultiIndex<INCREMENTER>::MultiIndex(bool first, ARGS&&... args)
    : m_incrementer(m_idx, first, std::forward<ARGS>(args)...)
    , m_first(first)
    , m_last(!first)
{}

template <class INCREMENTER>
std::size_t BRY::MultiIndex<INCREMENTER>::size() const {
    return m_idx.size();
}

template <class INCREMENTER>
INCREMENTER& BRY::MultiIndex<INCREMENTER>::incrementer() {
    return m_incrementer;
}

template <class INCREMENTER>
const INCREMENTER& BRY::MultiIndex<INCREMENTER>::incrementer() const {
    return m_incrementer;
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>& BRY::MultiIndex<INCREMENTER>::operator++() {
    if (m_incrementer.increment(m_idx)) {
        m_first = false;
        m_last = false;
    } else {
        m_last = true;
    }
    return *this;
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER> BRY::MultiIndex<INCREMENTER>::operator++(int) {
    BRY::MultiIndex<INCREMENTER> idx = *this;
    return ++idx;
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>& BRY::MultiIndex<INCREMENTER>::operator--() {
    if (m_incrementer.decrement(m_idx)) {
        m_first = false;
        m_last = false;
    } else {
        m_first = true;
    }
    return *this;
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER> BRY::MultiIndex<INCREMENTER>::operator--(int) {
    BRY::MultiIndex<INCREMENTER> idx = *this;
    return --idx;
}

template <class INCREMENTER>
bool BRY::MultiIndex<INCREMENTER>::first() const {
    return m_first;
}

template <class INCREMENTER>
bool BRY::MultiIndex<INCREMENTER>::last() const {
    return m_last;
}

template <class INCREMENTER>
BRY::bry_idx_t BRY::MultiIndex<INCREMENTER>::operator[](std::size_t d) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(d < m_idx.size(), "Subscript d is out of bounds");
#endif
    return m_idx[d];
}

template <class INCREMENTER>
std::vector<BRY::bry_idx_t>::const_iterator BRY::MultiIndex<INCREMENTER>::begin() const {
    return m_idx.begin();
}

template <class INCREMENTER>
std::vector<BRY::bry_idx_t>::const_iterator BRY::MultiIndex<INCREMENTER>::end() const {
    return m_idx.end();
}

template <class INCREMENTER>
const std::vector<BRY::bry_idx_t>& BRY::MultiIndex<INCREMENTER>::getIdxVector() const {
    return m_idx;
}

template <class INCREMENTER>
std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex<INCREMENTER>& idx) {
    os << BRY_LOG_BWHITE("[");
    for (std::size_t i = 0; i < idx.size(); ++i) {
        os << BRY_LOG_GREEN(idx[i]);
        if (i < idx.size() - 1)
            os << BRY_LOG_WHITE(", ");
    }
    os << BRY_LOG_BWHITE("]");
    return os;
}

BRY::MultiIndex<BRY::ExhaustiveIncrementer> BRY::mIdx(std::size_t sz, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementer>(true, sz, index_constraint);
}

BRY::MultiIndex<BRY::ExhaustiveIncrementer> BRY::rmIdx(std::size_t sz, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementer>(false, sz, index_constraint);
}

BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap> BRY::mIdxW(std::size_t sz, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap>(true, sz, index_constraint);
}

BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap> BRY::rmIdxW(std::size_t sz, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap>(false, sz, index_constraint);
}

BRY::MultiIndex<BRY::FixedNormIncrementer> BRY::mIdxFN(std::size_t sz, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::FixedNormIncrementer>(true, sz, index_constraint);
}

BRY::MultiIndex<BRY::FixedNormIncrementer> BRY::rmIdxFN(std::size_t sz, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::FixedNormIncrementer>(false, sz, index_constraint);
}

BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap> BRY::mIdxBEW(const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap>(true, index_bounds, index_constraint);
}

BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap> BRY::rmIdxBEW(const std::vector<bry_idx_t>& index_bounds, bry_idx_t index_constraint) {
    return BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap>(false, index_bounds, index_constraint);
}