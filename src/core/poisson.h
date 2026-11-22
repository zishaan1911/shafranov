#pragma once

#include "core/grid.h"

#include <vector>

namespace shafranov {

// Solves Delta* psi = rhs at interior grid points, with the boundary values of psi held fixed
// (they are read from psi on entry). rhs is only read at interior points. The discretisation is
// exactly delta_star's, so delta_star(psi) reproduces rhs at interior points to rounding error.
// Requires nZ - 1 to be a power of two (65, 129, ...): throws std::invalid_argument otherwise,
// or if psi or rhs has the wrong size.
void solve_poisson(const Grid& g, const std::vector<double>& rhs, std::vector<double>& psi);

} // namespace shafranov
