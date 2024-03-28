#pragma once

#include <cmath>
#include <algorithm>

#include "MultiIndex.h"

BRY::MultiIndex::MultiIndex(std::size_t sz, std::size_t l1_norm)
    : m_midx(sz * l1_norm, false)
    , m_l1_norm(l1_norm)
    , m_begin(true)
    , m_end(false)
{
    for (std::size_t i = 0; i < l1_norm; ++i) {
        m_midx[i] = true;
    }

    //for (bool b : m_midx) {
    //    std::cout << b << " ";
    //}
    //std::cout << std::endl;
}

std::size_t BRY::MultiIndex::size() const {
    return m_midx.size() / m_l1_norm;
}

BRY::MultiIndex& BRY::MultiIndex::operator++() {
    if (std::next_permutation(m_midx.begin(), m_midx.end(), std::greater())) {
        m_begin = false;
        m_end = false;
    } else {
        m_end = true;
        std::prev_permutation(m_midx.begin(), m_midx.end(), std::greater());
    }

    for (bool b : m_midx) {
        std::cout << b << " ";
    }
    std::cout << std::endl;
    return *this;
}

BRY::MultiIndex BRY::MultiIndex::operator++(int) {
    BRY::MultiIndex idx = *this;
    return ++idx;
}

bool BRY::MultiIndex::begin() const {
    return m_begin;
}

bool BRY::MultiIndex::end() const {
    return m_end;
}

std::size_t BRY::MultiIndex::operator[](std::size_t d) const {
#ifdef BRY_ENABLE_BOUNDS_CHECK
    ASSERT(d < size(), "Subscript d is out of bounds");
#endif
    std::size_t i_d = 0;
    for (std::size_t j = m_l1_norm * d; j < m_l1_norm * (d + 1); ++j) {
        i_d += m_midx[j];
    }
    return i_d;
}

std::ostream& operator<<(std::ostream& os, const BRY::MultiIndex& p) {
    std::string s = "[";
    os << BRY_LOG_BWHITE("[");
    for (std::size_t i = 0; i < p.size(); ++i) {
        os << BRY_LOG_GREEN(p[i]);
        if (i < p.size() - 1)
            os << BRY_LOG_WHITE(", ");
    }
    os << BRY_LOG_BWHITE("]");
    return os;
}