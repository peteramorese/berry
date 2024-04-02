#pragma once

#include "Bernstein.h"

template <std::size_t DIM>
Polynomial<DIM, BRY::Basis::Bernstein> BRY::BernsteinBasis<DIM>::to(const Polynomial<DIM, BRY::Basis::Power>& p, BRY::bry_deg_t degree_increase = 0) {
    
}
