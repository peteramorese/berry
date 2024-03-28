#pragma once

#include <cmath>
#include <algorithm>

#include "MultiIndex.h"

BRY::MultiIndex::MultiIndex(std::size_t sz, std::size_t l1_norm, bool begin)
    : m_combination(l1_norm + sz - 1, false)
    , m_idx(sz, 0)
    , m_l1_norm(l1_norm)
    , m_left(begin)
    , m_right(!begin)
{
    if (begin) {
        for (std::size_t i = 0; i < l1_norm; ++i) {
            m_combination[i] = true;
        }
        m_idx.front() = l1_norm;
    } else {
        auto it = m_combination.end();
        for (std::size_t i = 0; i < l1_norm; ++i) {
            *(--it) = true;
        }
        m_idx.back() = l1_norm;
    }
}

std::size_t BRY::MultiIndex::size() const {
    return m_idx.size();
}

std::size_t BRY::MultiIndex::l1Norm() const {
    return m_l1_norm;
}

BRY::MultiIndex& BRY::MultiIndex::operator++() {
    if (std::next_permutation(m_combination.begin(), m_combination.end(), std::greater())) {
        m_left = false;
        m_right = false;
    } else {
        m_right = true;
        std::prev_permutation(m_combination.begin(), m_combination.end(), std::greater());
    }

    updateIdx();
    return *this;
}

BRY::MultiIndex BRY::MultiIndex::operator++(int) {
    BRY::MultiIndex idx = *this;
    return ++idx;
}

BRY::MultiIndex& BRY::MultiIndex::operator--() {
    if (std::prev_permutation(m_combination.begin(), m_combination.end(), std::greater())) {
        m_left = false;
        m_right = false;
    } else {
        m_left = true;
        std::next_permutation(m_combination.begin(), m_combination.end(), std::greater());
    }

    updateIdx();
    return *this;
}

BRY::MultiIndex BRY::MultiIndex::operator--(int) {
    BRY::MultiIndex idx = *this;
    return --idx;
}

bool BRY::MultiIndex::left() const {
    return m_left;
}

bool BRY::MultiIndex::right() const {
    return m_right;
}

std::size_t BRY::MultiIndex::operator[](std::size_t d) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(d < m_idx.size(), "Subscript d is out of bounds");
#endif
    return m_idx[d];
}

std::vector<std::size_t>::const_iterator BRY::MultiIndex::begin() const {
    return m_idx.begin();
}

std::vector<std::size_t>::const_iterator BRY::MultiIndex::end() const {
    return m_idx.end();
}

void BRY::MultiIndex::updateIdx() {
    std::size_t d_i = 0;
    std::size_t sum = 0;
    for (bool bit : m_combination) {
        if (bit) {
            ++sum;
        } else {
            m_idx[d_i++] = sum;
            sum = 0;
        }
    }
    m_idx.back() = sum;
}

std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex& idx) {
    os << BRY_LOG_BWHITE("[");
    for (std::size_t i = 0; i < idx.size(); ++i) {
        os << BRY_LOG_GREEN(idx[i]);
        if (i < idx.size() - 1)
            os << BRY_LOG_WHITE(", ");
    }
    os << BRY_LOG_BWHITE("]");
    return os;
}