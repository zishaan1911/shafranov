#pragma once

#include <vector>

#include "core/grid.h"

namespace shafranov {

// Plasma current parametrisation: J_phi = lambda (beta0 R / R0 + (1 - beta0) R0 / R) shape(psi_n),
// shape = (1 - psi_n^alpha_m)^alpha_n, psi_n = 0 on the axis and 1 on the boundary.
struct Profile {
    double Ip;          // total toroidal plasma current (A); may be negative
    double beta0;       // in [0, 1]: share of the current carried by the pressure term
    double alpha_m = 1; // shape exponents (both > 0)
    double alpha_n = 2;
    double R0; // normalisation radius (m), > 0
};

// Fills jphi (resized to g.size()) with J_phi at grid points where inside[k] != 0, zero elsewhere,
// and returns lambda. Throws std::invalid_argument on size mismatches, psi_bndry == psi_axis, or
// invalid parameters, and std::runtime_error if the unscaled current integral is zero (no plasma
// points).
double plasma_current(const Grid& g, const std::vector<double>& psi, double psi_axis,
                      double psi_bndry, const std::vector<char>& inside, const Profile& p,
                      std::vector<double>& jphi);

// p'(psi_n) and FF'(psi_n) for a given lambda, as defined above (psi_n clamped to [0, 1]).
double pprime(const Profile& p, double lambda, double psi_n);
double ffprime(const Profile& p, double lambda, double psi_n);

} // namespace shafranov
