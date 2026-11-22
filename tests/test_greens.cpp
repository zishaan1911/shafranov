#include "core/greens.h"

#include "core/constants.h"

#include <doctest/doctest.h>

#include <array>
#include <cmath>

using shafranov::BField;
using shafranov::greens_b;
using shafranov::greens_psi;
using shafranov::mu0;
using shafranov::pi;

TEST_CASE("greens_psi is symmetric in source and field point") {
    const std::array<std::array<double, 4>, 4> pairs = {{
        {0.5, 0.0, 1.5, 0.3},
        {1.0, -0.4, 2.0, 0.7},
        {3.0, 1.0, 0.2, -1.0},
        {0.8, 0.1, 0.9, 0.2},
    }};
    for (const auto& p : pairs) {
        CHECK(greens_psi(p[0], p[1], p[2], p[3]) ==
              doctest::Approx(greens_psi(p[2], p[3], p[0], p[1])).epsilon(1e-12));
    }
}

TEST_CASE("on-axis field matches the analytic loop field") {
    const double Rc = 0.5;
    const double R = 1e-6;
    for (double d : {0.0, 0.3, 1.0}) {
        const BField b = greens_b(Rc, 0.0, R, d);
        const double expected = mu0 * Rc * Rc / (2.0 * std::pow(Rc * Rc + d * d, 1.5));
        CHECK(b.Z == doctest::Approx(expected).epsilon(1e-5));
        CHECK(std::abs(b.R) < 1e-9);
    }
}

TEST_CASE("far-field flux matches the magnetic dipole") {
    const double Rc = 0.5;
    const double m = pi * Rc * Rc;
    const double r = 200.0 * Rc;
    for (double theta : {pi / 3.0, pi / 6.0}) {
        const double R = r * std::sin(theta);
        const double Z = r * std::cos(theta);
        const double dipole = mu0 * m * R * R / (4.0 * pi * r * r * r);
        CHECK(greens_psi(Rc, 0.0, R, Z) == doctest::Approx(dipole).epsilon(1e-3));
    }
}

TEST_CASE("field equals derivatives of the flux") {
    const double Rc = 1.0;
    const double h = 1e-5;
    const std::array<std::array<double, 2>, 5> points = {{
        {0.5, 0.3},
        {1.5, -0.4},
        {2.0, 0.5},
        {0.8, -0.6},
        {1.2, 0.7},
    }};
    for (const auto& p : points) {
        const double R = p[0];
        const double Z = p[1];
        const double dpsi_dZ =
            (greens_psi(Rc, 0.0, R, Z + h) - greens_psi(Rc, 0.0, R, Z - h)) / (2.0 * h);
        const double dpsi_dR =
            (greens_psi(Rc, 0.0, R + h, Z) - greens_psi(Rc, 0.0, R - h, Z)) / (2.0 * h);
        const BField b = greens_b(Rc, 0.0, R, Z);
        CHECK(b.R == doctest::Approx(-dpsi_dZ / R).epsilon(1e-5));
        CHECK(b.Z == doctest::Approx(dpsi_dR / R).epsilon(1e-5));
    }
}

TEST_CASE("field is divergence free") {
    const double Rc = 1.0;
    const double h = 1e-4;
    const std::array<std::array<double, 2>, 3> points = {{
        {0.5, 0.3},
        {1.5, -0.4},
        {2.0, 0.5},
    }};
    for (const auto& p : points) {
        const double R = p[0];
        const double Z = p[1];
        const BField hi = greens_b(Rc, 0.0, R + h, Z);
        const BField lo = greens_b(Rc, 0.0, R - h, Z);
        const double d_rbr = ((R + h) * hi.R - (R - h) * lo.R) / (2.0 * h);
        const double dbz_dz = (greens_b(Rc, 0.0, R, Z + h).Z - greens_b(Rc, 0.0, R, Z - h).Z) /
                              (2.0 * h);
        const BField b = greens_b(Rc, 0.0, R, Z);
        const double b_mag = std::hypot(b.R, b.Z);
        CHECK(std::abs(d_rbr / R + dbz_dz) < 1e-6 * b_mag / h);
    }
}

TEST_CASE("mutual inductance of coaxial rings matches Maxwell's formula") {
    const double R1 = 1.0;
    const double R2 = 1.0;
    const double sep = 0.1;
    const double k = std::sqrt(4.0 * R1 * R2 / ((R1 + R2) * (R1 + R2) + sep * sep));
    const double K = std::comp_ellint_1(k);
    const double E = std::comp_ellint_2(k);
    const double expected = mu0 * std::sqrt(R1 * R2) * ((2.0 / k - k) * K - (2.0 / k) * E);
    CHECK(shafranov::mutual_inductance(R1, 0.0, R2, sep) ==
          doctest::Approx(expected).epsilon(1e-10));
}
