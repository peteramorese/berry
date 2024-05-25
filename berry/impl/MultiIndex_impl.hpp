#pragma once

#include "MultiIndex.h"

#include <cmath>
#include <algorithm>

/* Exhaustive Norm Incrementer */

BRY::ExhaustiveIncrementer::ExhaustiveIncrementer(bry_int_t*& _initial_idx, std::size_t sz, bool first, bry_int_t index_constraint)
    : m_linfty_norm(index_constraint)
{
    // If an external array is used, do not reallocate
    if (!_initial_idx) {
        _initial_idx = new bry_int_t[sz];
    }

    if (!first) {
        std::fill_n(_initial_idx, sz, index_constraint - 1);
    } else {
        std::fill_n(_initial_idx, sz, 0);
    }
}

BRY::bry_int_t BRY::ExhaustiveIncrementer::indexConstraint() const {
    return m_linfty_norm;
}


bool BRY::ExhaustiveIncrementer::increment(bry_int_t* current_idx, std::size_t sz) {
    ++current_idx[0];
    for (std::size_t i = 0; i < sz - 1; ++i) {
        if (current_idx[i] >= m_linfty_norm) {
            current_idx[i] = 0;
            ++current_idx[i + 1];
        }
    }
    return current_idx[sz - 1] < m_linfty_norm;
}

bool BRY::ExhaustiveIncrementer::decrement(bry_int_t* current_idx, std::size_t sz) {
    --current_idx[0];
    for (std::size_t i = 0; i < sz - 1; ++i) {
        if (current_idx[i] < 0) {
            current_idx[i] = m_linfty_norm - 1;
            --current_idx[i + 1];
        }
    }
    return current_idx[sz - 1] >= 0;
}

/* Exhaustive Norm Incrementer */

BRY::ExhaustiveIncrementerWrap::ExhaustiveIncrementerWrap(bry_int_t*& _initial_idx, std::size_t sz, bool first, bry_int_t index_constraint)
    : m_linfty_norm(index_constraint)
{
    // If an external array is used, do not reallocate
    if (!_initial_idx) {
        _initial_idx = new bry_int_t[sz];
    }

    if (!first) {
        std::fill_n(_initial_idx, sz, index_constraint - 1);
        m_wrapped_idx = std::pow(m_linfty_norm, sz) - 1;
    } else {
        std::fill_n(_initial_idx, sz, 0);
        m_wrapped_idx = 0;
    }
}

BRY::bry_int_t BRY::ExhaustiveIncrementerWrap::indexConstraint() const {
    return m_linfty_norm;
}


bool BRY::ExhaustiveIncrementerWrap::increment(bry_int_t* current_idx, std::size_t sz) {
    ++current_idx[0];
    ++m_wrapped_idx;
    for (std::size_t i = 0; i < sz - 1; ++i) {
        if (current_idx[i] >= m_linfty_norm) {
            current_idx[i] = 0;
            ++current_idx[i + 1];
        }
    }
    return current_idx[sz - 1] < m_linfty_norm;
}

bool BRY::ExhaustiveIncrementerWrap::decrement(bry_int_t* current_idx, std::size_t sz) {
    --current_idx[0];
    --m_wrapped_idx;
    for (std::size_t i = 0; i < sz - 1; ++i) {
        if (current_idx[i] < 0) {
            current_idx[i] = m_linfty_norm - 1;
            --current_idx[i + 1];
        }
    }
    return current_idx[sz - 1] >= 0;
}

BRY::bry_int_t BRY::ExhaustiveIncrementerWrap::wrappedIdx() const {
    return m_wrapped_idx;
}

/* Fixed Norm Incrementer */

BRY::FixedNormIncrementer::FixedNormIncrementer(bry_int_t*& _initial_idx, std::size_t sz, bool first, bry_int_t index_constraint)
    : m_combination(index_constraint + sz - 2, false)
    , m_l1_norm(index_constraint)
{
    // If an external array is used, do not reallocate
    if (!_initial_idx) {
        _initial_idx = new bry_int_t[sz];
    }

    if (first) {
        std::fill_n(m_combination.begin(), index_constraint - 1, true);
        _initial_idx[0] = index_constraint - 1;
    } else {
        auto it = m_combination.end();
        for (std::size_t i = 0; i < index_constraint - 1; ++i) {
            *(--it) = true;
        }
        _initial_idx[sz - 1] = index_constraint - 1;
    }
}

BRY::bry_int_t BRY::FixedNormIncrementer::indexConstraint() const {
    return m_l1_norm;
}

bool BRY::FixedNormIncrementer::increment(bry_int_t* current_idx, std::size_t sz) {
    if (std::next_permutation(m_combination.begin(), m_combination.end(), std::greater())) {
        updateIdx(current_idx, sz);
        return true;
    } else {
        std::prev_permutation(m_combination.begin(), m_combination.end(), std::greater());
        return false;
    }
}

bool BRY::FixedNormIncrementer::decrement(bry_int_t* current_idx, std::size_t sz) {
    if (std::prev_permutation(m_combination.begin(), m_combination.end(), std::greater())) {
        updateIdx(current_idx, sz);
        return true;
    } else {
        std::next_permutation(m_combination.begin(), m_combination.end(), std::greater());
        return false;
    }
}

void BRY::FixedNormIncrementer::updateIdx(bry_int_t* current_idx, std::size_t sz) const {
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
    current_idx[sz - 1] = sum;
}

/* Bounded Exhaustive Incrementer Wrap */

BRY::BoundedExhaustiveIncrementerWrap::BoundedExhaustiveIncrementerWrap(bry_int_t*& _initial_idx, std::size_t sz, bool first, const std::vector<bry_int_t>& index_bounds, bry_int_t index_constraint) 
    : m_index_bounds(index_bounds)
    , m_linfty_norm(index_constraint)
{
    // If an external array is used, do not reallocate
    if (!_initial_idx) {
        _initial_idx = new bry_int_t[index_bounds.size()];
    }

    if (first) {
        m_wrapped_idx = 0;
        std::fill_n(_initial_idx, sz, 0);
    } else {
        m_wrapped_idx = 0;
        for (std::size_t i = 0; i < sz; ++i) {
            _initial_idx[i] = index_bounds[i] - 1;
            m_wrapped_idx += std::pow(m_linfty_norm, i) * (index_bounds[i] - 1);
        }
    }
}

const std::vector<BRY::bry_int_t>& BRY::BoundedExhaustiveIncrementerWrap::indexConstraint() const {
    return m_index_bounds;
}

bool BRY::BoundedExhaustiveIncrementerWrap::increment(bry_int_t* current_idx, std::size_t sz) {
    ++current_idx[0];
    ++m_wrapped_idx;
    for (std::size_t i = 0; i < sz - 1; ++i) {
        if (current_idx[i] >= m_index_bounds[i]) {
            current_idx[i] = 0;
            ++current_idx[i + 1];
            m_wrapped_idx += std::pow(m_linfty_norm, i) * (m_linfty_norm - m_index_bounds[i]);
        }
    }
    return current_idx[sz - 1] < m_index_bounds.back();
}

bool BRY::BoundedExhaustiveIncrementerWrap::decrement(bry_int_t* current_idx, std::size_t sz) {
    --current_idx[0];
    --m_wrapped_idx;
    for (std::size_t i = 0; i < sz - 1; ++i) {
        if (current_idx[i] < 0) {
            current_idx[i] = m_index_bounds[i] - 1;
            --current_idx[i + 1];
            m_wrapped_idx -= std::pow(m_linfty_norm, i) * (m_linfty_norm - m_index_bounds[i]);
        }
    }
    return current_idx[sz - 1] >= 0;
}

BRY::bry_int_t BRY::BoundedExhaustiveIncrementerWrap::wrappedIdx() const {
    return m_wrapped_idx;
}

/* MultiIndex */

template <class INCREMENTER>
template <typename ... ARGS>
BRY::MultiIndex<INCREMENTER>::MultiIndex(std::size_t sz, bool first, ARGS&&... args)
    : m_idx(nullptr)
    , m_sz(sz)
    , m_incrementer(m_idx, m_sz, first, std::forward<ARGS>(args)...)
    , m_first(first)
    , m_last(!first)
    , m_external_arr(false)
{}

template <class INCREMENTER>
template <typename ... ARGS>
BRY::MultiIndex<INCREMENTER>::MultiIndex(bry_int_t* array, std::size_t sz, bool first, ARGS&&... args)
    : m_idx(array)
    , m_sz(sz)
    , m_incrementer(m_idx, m_sz, first, std::forward<ARGS>(args)...)
    , m_first(first)
    , m_last(!first)
    , m_external_arr(true)
{
    ASSERT(sz, "Size must be geq 0");
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>::~MultiIndex() {
    if (!m_external_arr)
        delete[] m_idx;
}

template <class INCREMENTER>
std::size_t BRY::MultiIndex<INCREMENTER>::size() const {
    return m_sz;
}

template <class INCREMENTER>
INCREMENTER& BRY::MultiIndex<INCREMENTER>::inc() {
    return m_incrementer;
}

template <class INCREMENTER>
const INCREMENTER& BRY::MultiIndex<INCREMENTER>::inc() const {
    return m_incrementer;
}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>& BRY::MultiIndex<INCREMENTER>::operator++() {
    if (m_incrementer.increment(m_idx, m_sz)) {
        m_first = false;
        m_last = false;
    } else {
        m_last = true;
    }
    return *this;
}

//template <class INCREMENTER>
//BRY::MultiIndex<INCREMENTER> BRY::MultiIndex<INCREMENTER>::operator++(int) {
//    BRY::MultiIndex<INCREMENTER> idx = *this;
//    return ++idx;
//}

template <class INCREMENTER>
BRY::MultiIndex<INCREMENTER>& BRY::MultiIndex<INCREMENTER>::operator--() {
    if (m_incrementer.decrement(m_idx, m_sz)) {
        m_first = false;
        m_last = false;
    } else {
        m_first = true;
    }
    return *this;
}

//template <class INCREMENTER>
//BRY::MultiIndex<INCREMENTER> BRY::MultiIndex<INCREMENTER>::operator--(int) {
//    BRY::MultiIndex<INCREMENTER> idx = *this;
//    return --idx;
//}

template <class INCREMENTER>
bool BRY::MultiIndex<INCREMENTER>::first() const {
    return m_first;
}

template <class INCREMENTER>
bool BRY::MultiIndex<INCREMENTER>::last() const {
    return m_last;
}

template <class INCREMENTER>
BRY::bry_int_t BRY::MultiIndex<INCREMENTER>::operator[](std::size_t d) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(d < m_sz, "Subscript d is out of bounds");
#endif
    return m_idx[d];
}

template <class INCREMENTER>
const BRY::bry_int_t* BRY::MultiIndex<INCREMENTER>::begin() const {
    return m_idx;
}

template <class INCREMENTER>
const BRY::bry_int_t* BRY::MultiIndex<INCREMENTER>::end() const {
    return m_idx + m_sz;
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

BRY::MultiIndex<BRY::ExhaustiveIncrementer> BRY::mIdx(std::size_t sz, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementer>(sz, true, index_constraint);
}

BRY::MultiIndex<BRY::ExhaustiveIncrementer> BRY::rmIdx(std::size_t sz, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementer>(sz, false, index_constraint);
}

BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap> BRY::mIdxW(std::size_t sz, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap>(sz, true, index_constraint);
}

BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap> BRY::rmIdxW(std::size_t sz, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::ExhaustiveIncrementerWrap>(sz, false, index_constraint);
}

BRY::MultiIndex<BRY::FixedNormIncrementer> BRY::mIdxFN(std::size_t sz, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::FixedNormIncrementer>(sz, true, index_constraint);
}

BRY::MultiIndex<BRY::FixedNormIncrementer> BRY::rmIdxFN(std::size_t sz, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::FixedNormIncrementer>(sz, false, index_constraint);
}

BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap> BRY::mIdxBEW(const std::vector<bry_int_t>& index_bounds, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap>(index_bounds.size(), true, index_bounds, index_constraint);
}

BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap> BRY::rmIdxBEW(const std::vector<bry_int_t>& index_bounds, bry_int_t index_constraint) {
    return BRY::MultiIndex<BRY::BoundedExhaustiveIncrementerWrap>(index_bounds.size(), false, index_bounds, index_constraint);
}