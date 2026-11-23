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

// Complete elliptic integrals K and E of modulus sqrt(k2), together, by the arithmetic-geometric
// mean: several times faster than std::comp_ellint_1 and _2, and as accurate.
void ellipke(double k2, double& K, double& E) {
    double a = 1.0;
    double b = std::sqrt(1.0 - k2);
    double sum = 0.5 * k2; // sum of 2^(n-1) c_n^2, with c_0^2 = k2
    double weight = 0.5;
    for (int n = 0; n < 32; ++n) {
        const double c = 0.5 * (a - b);
        const double next = 0.5 * (a + b);
        b = std::sqrt(a * b);
        a = next;
        weight *= 2.0;
        sum += weight * c * c;
        if (std::abs(c) <= 1e-16 * a) {
            break;
        }
    }
    K = pi / (2.0 * a);
    E = K * (1.0 - sum);
}

} // namespace

double greens_psi(double Rc, double Zc, double R, double Z) {
    const double k2 = clamped_k2(Rc, R, Z - Zc);
    const double k = std::sqrt(k2);
    double K, E;
    ellipke(k2, K, E);
    return (mu0 / (2.0 * pi)) * std::sqrt(R * Rc) * ((2.0 - k2) * K - 2.0 * E) / k;
}

BField greens_b(double Rc, double Zc, double R, double Z) {
    const double d = Z - Zc;
    const double k2 = clamped_k2(Rc, R, d);
    double K, E;
    ellipke(k2, K, E);
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
