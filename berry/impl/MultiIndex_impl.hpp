#pragma once

#include "MultiIndex.h"

#include <cmath>
#include <algorithm>

/* Fixed Norm Incrementer */

BRY::FixedNormIncrementer::FixedNormIncrementer(std::size_t sz, bry_idx_t index_constraint, bool begin, std::vector<bry_idx_t>& initial_idx)
    : m_combination(index_constraint + sz - 2, false)
    , m_l1_norm(index_constraint)
{
    if (begin) {
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

/* Exhaustive Norm Incrementer */

BRY::ExhaustiveIncrementer::ExhaustiveIncrementer(std::size_t sz, bry_idx_t index_constraint, bool begin, std::vector<bry_idx_t>& initial_idx)
    : m_linfty_norm(index_constraint)
{
    if (!begin) {
        for (bry_idx_t& i : initial_idx)
            i = index_constraint - 1;
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

/* MultiIndex */

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>::MultiIndex(std::size_t sz, bry_idx_t index_constraint, bool begin)
    : m_idx(sz, 0)
    , m_incrementer(sz, index_constraint, begin, m_idx)
    , m_left(begin)
    , m_right(!begin)
{}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>::MultiIndex(INCREMENTER&& incrementer, std::vector<bry_idx_t>&& initial_idx) 
    : m_idx(std::move(initial_idx))
    , m_incrementer(std::move(incrementer))
    , m_left(begin)
    , m_right(!begin)
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
        m_left = false;
        m_right = false;
    } else {
        m_right = true;
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
        m_left = false;
        m_right = false;
    } else {
        m_left = true;
    }
    return *this;
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER> BRY::MultiIndex<INCREMENTER>::operator--(int) {
    BRY::MultiIndex<INCREMENTER> idx = *this;
    return --idx;
}

template <class INCREMENTER>
bool BRY::MultiIndex<INCREMENTER>::left() const {
    return m_left;
}

template <class INCREMENTER>
bool BRY::MultiIndex<INCREMENTER>::right() const {
    return m_right;
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