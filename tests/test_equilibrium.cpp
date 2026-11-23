#include "core/constants.h"
#include "core/equilibrium.h"
#include "core/greens.h"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace shafranov;

namespace {

constexpr double Ip = 1e5;

std::vector<Point> box(double Rlo, double Rhi, double Zlo, double Zhi) {
    const int n = 40;
    std::vector<Point> pts;
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rlo + (Rhi - Rlo) * k / n, Zlo});
    }
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rhi, Zlo + (Zhi - Zlo) * k / n});
    }
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rhi - (Rhi - Rlo) * k / n, Zhi});
    }
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rlo, Zhi - (Zhi - Zlo) * k / n});
    }
    return pts;
}

Machine test_machine() {
    Machine m;
    m.coils.push_back({"up", "V", 2.0, 1.0, 0.05, 0.05, 1.0, 0.0});
    m.coils.push_back({"down", "V", 2.0, -1.0, 0.05, 0.05, 1.0, 0.0});
    m.limiter = box(0.6, 1.2, -0.4, 0.4);
    return m;
}

std::size_t at(const Grid& g, int i, int j) { return static_cast<std::size_t>(g.idx(i, j)); }

// Coil current that balances the plasma's outward hoop force at (0.9, 0) for Ip = 1e5 A.
double balance_current() {
    const double R0 = 0.9;
    const double a = 0.2;
    const double bz_req = -mu0 * Ip / (4.0 * pi * R0) * (std::log(8.0 * R0 / a) + 0.3 + 0.5 - 1.5);
    const double bz_per_amp = greens_b(2.0, 1.0, R0, 0.0).Z + greens_b(2.0, -1.0, R0, 0.0).Z;
    return bz_req / bz_per_amp;
}

struct Fixture {
    Grid grid{65, 65, 0.5, 1.3, -0.5, 0.5};
    MachineGrid mg = make_machine_grid(test_machine(), grid);
    Profile profile{Ip, 0.3, 1.0, 2.0, 0.9};
    InitialPlasma initial{0.9, 0.0, 0.2};
    SolveOptions options{};
};

} // namespace

TEST_CASE("equilibrium converges to a centred, consistent solution") {
    Fixture f;
    Currents currents{{balance_current()}, {}};
    const Equilibrium eq = solve_equilibrium(f.mg, currents, f.profile, f.initial, f.options);
    const Grid& g = f.grid;

    // 1. convergence and axis position
    CHECK(eq.converged);
    CHECK(std::abs(eq.topology.axis.R - 0.9) < 0.15);
    CHECK(std::abs(eq.topology.axis.Z) < 1e-6);

    // 2. total current
    double total = 0.0;
    double max_rhs = 0.0;
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            total += eq.jphi[at(g, i, j)] * g.dR() * g.dZ();
            max_rhs = std::max(max_rhs, std::abs(mu0 * g.R(i) * eq.jphi[at(g, i, j)]));
        }
    }
    CHECK(total == doctest::Approx(Ip).epsilon(1e-9));

    // 3. Grad-Shafranov consistency of the plasma flux
    std::vector<double> ds;
    delta_star(g, eq.psi_plasma, ds);
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            const double expected = -mu0 * g.R(i) * eq.jphi[at(g, i, j)];
            CHECK(std::abs(ds[at(g, i, j)] - expected) < 1e-6 * max_rhs);
        }
    }

    // 4. free-space flux at nodes well outside the plasma
    std::vector<std::pair<int, int>> far;
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            bool clear = true;
            for (int di = -3; di <= 3 && clear; ++di) {
                for (int dj = -3; dj <= 3 && clear; ++dj) {
                    const int ii = i + di;
                    const int jj = j + dj;
                    if (ii >= 0 && ii < g.nR && jj >= 0 && jj < g.nZ &&
                        eq.jphi[at(g, ii, jj)] != 0.0) {
                        clear = false;
                    }
                }
            }
            if (clear) {
                far.emplace_back(i, j);
            }
        }
    }
    REQUIRE(far.size() >= 3);
    const std::vector<std::size_t> picks = {0, far.size() / 2, far.size() - 1};
    for (const std::size_t p : picks) {
        const auto [i, j] = far[p];
        double direct = 0.0;
        for (int a = 0; a < g.nR; ++a) {
            for (int b = 0; b < g.nZ; ++b) {
                const double J = eq.jphi[at(g, a, b)];
                if (J != 0.0) {
                    direct += greens_psi(g.R(a), g.Z(b), g.R(i), g.Z(j)) * J * g.dR() * g.dZ();
                }
            }
        }
        CHECK(eq.psi_plasma[at(g, i, j)] == doctest::Approx(direct).epsilon(1e-2));
    }
}

TEST_CASE("reversing the plasma and coil currents gives the same axis") {
    Fixture f;
    Currents currents{{balance_current()}, {}};
    const Equilibrium eq = solve_equilibrium(f.mg, currents, f.profile, f.initial, f.options);

    Profile reversed = f.profile;
    reversed.Ip = -Ip;
    Currents flipped{{-balance_current()}, {}};
    const Equilibrium eq2 = solve_equilibrium(f.mg, flipped, reversed, f.initial, f.options);

    CHECK(eq2.converged);
    CHECK(eq2.topology.axis.R == doctest::Approx(eq.topology.axis.R).epsilon(1e-6));
    CHECK(std::abs(eq2.topology.axis.Z - eq.topology.axis.Z) < 1e-6);
}

TEST_CASE("feedback is called once per iteration") {
    Fixture f;
    Currents currents{{balance_current()}, {}};
    int calls = 0;
    const Feedback count = [&calls](const Equilibrium&, Currents&) { ++calls; };
    const Equilibrium eq =
        solve_equilibrium(f.mg, currents, f.profile, f.initial, f.options, count);
    CHECK(eq.iterations > 0);
    CHECK(calls == eq.iterations);
}

TEST_CASE("solve_equilibrium rejects a wrong number of circuit currents") {
    Fixture f;
    Currents currents{{1.0, 2.0}, {}};
    CHECK_THROWS_AS(solve_equilibrium(f.mg, currents, f.profile, f.initial, f.options),
                    std::invalid_argument);
}
