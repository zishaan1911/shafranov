#include "core/greens.h"

#include "core/constants.h"

#include <algorithm>
#include <cmath>

namespace shafranov {

namespace {

// Keeps k away from 0 and 1 so evaluating on the filament itself stays finite.
double clamped_k2(double Rc, double R, double d) {
    const double k2 = 4.0 * R * Rc / ((R + Rc) * (R + Rc) + d * d);
    return std::clamp(k2, 1e-10, 1.0 - 1e-10);
}

} // namespace

double greens_psi(double Rc, double Zc, double R, double Z) {
    const double k2 = clamped_k2(Rc, R, Z - Zc);
    const double k = std::sqrt(k2);
    const double K = std::comp_ellint_1(k);
    const double E = std::comp_ellint_2(k);
    return (mu0 / (2.0 * pi)) * std::sqrt(R * Rc) * ((2.0 - k2) * K - 2.0 * E) / k;
}

BField greens_b(double Rc, double Zc, double R, double Z) {
    const double d = Z - Zc;
    const double k2 = clamped_k2(Rc, R, d);
    const double k = std::sqrt(k2);
    const double K = std::comp_ellint_1(k);
    const double E = std::comp_ellint_2(k);
    const double alpha2 = (R - Rc) * (R - Rc) + d * d;
    const double beta = std::sqrt((R + Rc) * (R + Rc) + d * d);
    const double c = mu0 / (2.0 * pi);
    return BField{
        c * d / (R * beta) * (-K + (Rc * Rc + R * R + d * d) / alpha2 * E),
        c / beta * (K + (Rc * Rc - R * R - d * d) / alpha2 * E),
    };
}

double mutual_inductance(double R1, double Z1, double R2, double Z2) {
    return 2.0 * pi * greens_psi(R1, Z1, R2, Z2);
}

} // namespace shafranov
