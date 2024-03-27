#pragma once

#include <cmath>
//#include <experimental/array>

#include "MultiIndex.h"

BRY::MultiIndex::MultiIndex(std::size_t sz, std::size_t l1_norm)
    : m_midx(sz)
    , m_l1_norm(l1_norm)
    , m_i(0ul)
{
    m_midx[0] = l1_norm;
}

BRY::MultiIndex& BRY::MultiIndex::operator++() {

}

BRY::MultiIndex BRY::MultiIndex::operator++(int) {

}

bool BRY::MultiIndex::end() const {

}

std::size_t BRY::MultiIndex::operator[](std::size_t i) const {
    return m_midx[i];
}