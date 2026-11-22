#include "core/grid.h"
#include "core/poisson.h"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace {

using shafranov::Grid;

std::size_t at(const Grid& g, int i, int j) { return static_cast<std::size_t>(g.idx(i, j)); }

bool is_boundary(const Grid& g, int i, int j) {
    return i == 0 || j == 0 || i == g.nR - 1 || j == g.nZ - 1;
}

// Smooth rhs everywhere; smooth boundary values, zero start in the interior.
void smooth_problem(const Grid& g, std::vector<double>& rhs, std::vector<double>& psi) {
    rhs.assign(static_cast<std::size_t>(g.size()), 0.0);
    psi.assign(static_cast<std::size_t>(g.size()), 0.0);
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            const double R = g.R(i);
            const double Z = g.Z(j);
            rhs[at(g, i, j)] = std::sin(3.0 * R + 1.0) * std::cos(2.0 * Z) + 0.5 * R * Z + 2.0;
            if (is_boundary(g, i, j)) {
                psi[at(g, i, j)] = std::cos(2.0 * R - Z) + 0.3 * R * R * Z + 1.0;
            }
        }
    }
}

void check_discrete_inverse(const Grid& g) {
    std::vector<double> rhs;
    std::vector<double> psi;
    smooth_problem(g, rhs, psi);
    const std::vector<double> boundary = psi;

    shafranov::solve_poisson(g, rhs, psi);

    double max_rhs = 0.0;
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            max_rhs = std::max(max_rhs, std::abs(rhs[at(g, i, j)]));
        }
    }

    std::vector<double> out;
    shafranov::delta_star(g, psi, out);
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            if (is_boundary(g, i, j)) {
                CHECK(psi[at(g, i, j)] == boundary[at(g, i, j)]);
            } else {
                CHECK(std::abs(out[at(g, i, j)] - rhs[at(g, i, j)]) < 1e-8 * max_rhs);
            }
        }
    }
}

double solovev_error(int n, double* max_psi) {
    const double a = 1.3;
    const double b = -0.7;
    const double c = 0.4;
    const Grid g(n, n, 0.5, 1.5, -1.0, 1.0);
    std::vector<double> exact(static_cast<std::size_t>(g.size()));
    std::vector<double> rhs(static_cast<std::size_t>(g.size()));
    std::vector<double> psi(static_cast<std::size_t>(g.size()), 0.0);
    double peak = 0.0;
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            const double R = g.R(i);
            const double Z = g.Z(j);
            const double p = a * R * R * R * R / 8.0 + b * Z * Z / 2.0 + c * R * R;
            exact[at(g, i, j)] = p;
            rhs[at(g, i, j)] = a * R * R + b;
            peak = std::max(peak, std::abs(p));
            if (is_boundary(g, i, j)) {
                psi[at(g, i, j)] = p;
            }
        }
    }
    shafranov::solve_poisson(g, rhs, psi);
    double err = 0.0;
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            err = std::max(err, std::abs(psi[at(g, i, j)] - exact[at(g, i, j)]));
        }
    }
    if (max_psi != nullptr) {
        *max_psi = peak;
    }
    return err;
}

double gaussian_error(int n) {
    const double s = 0.1;
    const Grid g(n, n, 0.5, 1.5, -1.0, 1.0);
    std::vector<double> exact(static_cast<std::size_t>(g.size()));
    std::vector<double> rhs(static_cast<std::size_t>(g.size()));
    std::vector<double> psi(static_cast<std::size_t>(g.size()), 0.0);
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            const double R = g.R(i);
            const double Z = g.Z(j);
            const double dr = R - 1.0;
            const double f = std::exp(-(dr * dr + Z * Z) / s);
            // f_R = -2 dr f / s, f_RR = (4 dr^2 / s^2 - 2 / s) f, likewise in Z.
            const double f_R = -2.0 * dr / s * f;
            const double f_RR = (4.0 * dr * dr / (s * s) - 2.0 / s) * f;
            const double f_ZZ = (4.0 * Z * Z / (s * s) - 2.0 / s) * f;
            exact[at(g, i, j)] = f;
            rhs[at(g, i, j)] = f_RR - f_R / R + f_ZZ;
            if (is_boundary(g, i, j)) {
                psi[at(g, i, j)] = f;
            }
        }
    }
    shafranov::solve_poisson(g, rhs, psi);
    double err = 0.0;
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            err = std::max(err, std::abs(psi[at(g, i, j)] - exact[at(g, i, j)]));
        }
    }
    return err;
}

} // namespace

TEST_CASE("solve_poisson inverts delta_star exactly on a 33 x 33 grid") {
    check_discrete_inverse(Grid(33, 33, 0.5, 1.5, -1.0, 1.0));
}

TEST_CASE("solve_poisson recovers a Solov'ev equilibrium at second order") {
    double peak = 0.0;
    const double e33 = solovev_error(33, nullptr);
    const double e65 = solovev_error(65, &peak);
    CHECK(e33 / e65 > 3.5);
    CHECK(e65 < 1e-4 * peak);
}

TEST_CASE("solve_poisson converges at second order for a smooth field") {
    const double e33 = gaussian_error(33);
    const double e65 = gaussian_error(65);
    const double e129 = gaussian_error(129);
    CHECK(e33 / e65 > 3.5);
    CHECK(e65 / e129 > 3.5);
}

TEST_CASE("solve_poisson works on non-square grids") {
    check_discrete_inverse(Grid(40, 65, 0.5, 1.5, -1.0, 1.0));
}

TEST_CASE("solve_poisson rejects bad sizes") {
    const Grid bad(33, 64, 0.5, 1.5, -1.0, 1.0);
    std::vector<double> rhs(static_cast<std::size_t>(bad.size()), 0.0);
    std::vector<double> psi(static_cast<std::size_t>(bad.size()), 0.0);
    CHECK_THROWS_AS(shafranov::solve_poisson(bad, rhs, psi), std::invalid_argument);

    const Grid g(33, 33, 0.5, 1.5, -1.0, 1.0);
    std::vector<double> psi_ok(static_cast<std::size_t>(g.size()), 0.0);
    CHECK_THROWS_AS(shafranov::solve_poisson(g, std::vector<double>(5), psi_ok),
                    std::invalid_argument);
}
