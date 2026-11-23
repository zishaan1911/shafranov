#include "core/equilibrium.h"

#include "core/constants.h"
#include "core/greens.h"
#include "core/poisson.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <thread>

namespace shafranov {

namespace {

// Runs f(k) for k in [0, count) on every hardware thread; each k is independent.
template <class F> void parallel_for(std::size_t count, F f) {
    const std::size_t workers = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    for (std::size_t w = 0; w < workers; ++w) {
        threads.emplace_back([=] {
            for (std::size_t k = w; k < count; k += workers) {
                f(k);
            }
        });
    }
    for (std::thread& t : threads) {
        t.join();
    }
}

} // namespace

MachineGrid make_machine_grid(const Machine& m, const Grid& g) {
    MachineGrid mg{g, circuit_names(m), circuit_windings(m), {}, {}, {}, {}, m.limiter, m.vessel};
    const auto n = static_cast<std::size_t>(g.size());

    mg.circuit_psi.assign(mg.circuits.size(), std::vector<double>(n, 0.0));
    mg.vessel_psi.assign(m.vessel.size(), std::vector<double>(n, 0.0));
    parallel_for(n, [&](std::size_t k) {
        const double R = g.R(static_cast<int>(k) / g.nZ);
        const double Z = g.Z(static_cast<int>(k) % g.nZ);
        for (std::size_t c = 0; c < mg.circuits.size(); ++c) {
            for (const Winding& w : mg.windings[c]) {
                mg.circuit_psi[c][k] += w.turns * greens_psi(w.R, w.Z, R, Z);
            }
        }
        for (std::size_t v = 0; v < m.vessel.size(); ++v) {
            mg.vessel_psi[v][k] = greens_psi(m.vessel[v].R, m.vessel[v].Z, R, Z);
        }
    });

    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            if (i == 0 || i == g.nR - 1 || j == 0 || j == g.nZ - 1) {
                mg.boundary.push_back(g.idx(i, j));
            }
        }
    }

    mg.boundary_greens.resize(mg.boundary.size() * n);
    parallel_for(mg.boundary.size(), [&](std::size_t b) {
        const double Rb = g.R(mg.boundary[b] / g.nZ);
        const double Zb = g.Z(mg.boundary[b] % g.nZ);
        for (int i = 0; i < g.nR; ++i) {
            for (int j = 0; j < g.nZ; ++j) {
                mg.boundary_greens[b * n + static_cast<std::size_t>(g.idx(i, j))] =
                    greens_psi(g.R(i), g.Z(j), Rb, Zb);
            }
        }
    });
    return mg;
}

Equilibrium solve_equilibrium(const MachineGrid& mg, Currents& currents, const Profile& profile,
                              const InitialPlasma& initial, const SolveOptions& options,
                              const Feedback& feedback) {
    if (currents.circuits.size() != mg.circuits.size()) {
        throw std::invalid_argument("solve_equilibrium: wrong number of circuit currents");
    }
    if (!currents.vessel.empty() && currents.vessel.size() != mg.vessel_psi.size()) {
        throw std::invalid_argument("solve_equilibrium: wrong number of vessel currents");
    }
    if (options.max_iterations <= 0 || !(options.relax > 0.0 && options.relax <= 1.0) ||
        !(options.tolerance > 0.0)) {
        throw std::invalid_argument("solve_equilibrium: invalid options");
    }
    if (!(initial.a > 0.0)) {
        throw std::invalid_argument("solve_equilibrium: initial plasma radius must be positive");
    }

    const Grid& g = mg.grid;
    const auto n = static_cast<std::size_t>(g.size());
    const double cell = g.dR() * g.dZ();

    auto external = [&](const Currents& c) {
        std::vector<double> psi(n, 0.0);
        for (std::size_t k = 0; k < c.circuits.size(); ++k) {
            for (std::size_t p = 0; p < n; ++p) {
                psi[p] += c.circuits[k] * mg.circuit_psi[k][p];
            }
        }
        for (std::size_t k = 0; k < c.vessel.size(); ++k) {
            for (std::size_t p = 0; p < n; ++p) {
                psi[p] += c.vessel[k] * mg.vessel_psi[k][p];
            }
        }
        return psi;
    };

    auto plasma_flux = [&](const std::vector<double>& jphi) {
        std::vector<double> rhs(n);
        for (int i = 0; i < g.nR; ++i) {
            for (int j = 0; j < g.nZ; ++j) {
                const auto k = static_cast<std::size_t>(g.idx(i, j));
                rhs[k] = -mu0 * g.R(i) * jphi[k];
            }
        }
        std::vector<std::size_t> plasma_nodes;
        for (std::size_t k = 0; k < n; ++k) {
            if (jphi[k] != 0.0) {
                plasma_nodes.push_back(k);
            }
        }
        std::vector<double> psi(n, 0.0);
        for (std::size_t b = 0; b < mg.boundary.size(); ++b) {
            const double* row = &mg.boundary_greens[b * n];
            double sum = 0.0;
            for (const std::size_t k : plasma_nodes) {
                sum += row[k] * jphi[k];
            }
            psi[static_cast<std::size_t>(mg.boundary[b])] = sum * cell;
        }
        solve_poisson(g, rhs, psi);
        return psi;
    };

    std::vector<double> jphi(n, 0.0);
    double total = 0.0;
    if (!initial.jphi.empty()) {
        if (initial.jphi.size() != n) {
            throw std::invalid_argument("solve_equilibrium: initial jphi has the wrong size");
        }
        jphi = initial.jphi;
        for (const double J : jphi) {
            total += J;
        }
    } else {
        for (int i = 0; i < g.nR; ++i) {
            for (int j = 0; j < g.nZ; ++j) {
                const double dr = g.R(i) - initial.R;
                const double dz = g.Z(j) - initial.Z;
                const double r2 = (dr * dr + dz * dz) / (initial.a * initial.a);
                if (r2 < 1.0) {
                    jphi[static_cast<std::size_t>(g.idx(i, j))] = 1.0 - r2;
                    total += 1.0 - r2;
                }
            }
        }
    }
    if (total == 0.0) {
        throw std::invalid_argument("solve_equilibrium: initial plasma carries no current");
    }
    for (double& J : jphi) {
        J *= profile.Ip / (total * cell);
    }
    std::vector<double> psi_pl = plasma_flux(jphi);

    Equilibrium eq;
    for (int it = 1; it <= options.max_iterations; ++it) {
        auto total_psi = [&]() {
            const std::vector<double> ext = external(currents);
            std::vector<double> psi(n);
            for (std::size_t k = 0; k < n; ++k) {
                psi[k] = psi_pl[k] + ext[k];
            }
            return psi;
        };

        eq.psi = total_psi();
        if (feedback) {
            eq.psi_plasma = psi_pl;
            eq.jphi = jphi;
            feedback(eq, currents);
            if (currents.circuits.size() != mg.circuits.size() ||
                (!currents.vessel.empty() && currents.vessel.size() != mg.vessel_psi.size())) {
                throw std::invalid_argument("solve_equilibrium: feedback resized the currents");
            }
            eq.psi = total_psi();
        }
        eq.topology = find_topology(g, eq.psi, mg.limiter);
        eq.lambda = plasma_current(g, eq.psi, eq.topology.axis.psi, eq.topology.psi_bndry,
                                   eq.topology.inside, profile, jphi);

        const std::vector<double> new_pl = plasma_flux(jphi);
        double change = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            change = std::max(change, std::abs(new_pl[k] - psi_pl[k]));
        }
        change /= std::abs(eq.topology.axis.psi - eq.topology.psi_bndry);
        for (std::size_t k = 0; k < n; ++k) {
            psi_pl[k] = options.relax * new_pl[k] + (1.0 - options.relax) * psi_pl[k];
        }
        eq.iterations = it;

        if (change < options.tolerance) {
            eq.psi = total_psi();
            eq.topology = find_topology(g, eq.psi, mg.limiter);
            eq.lambda = plasma_current(g, eq.psi, eq.topology.axis.psi, eq.topology.psi_bndry,
                                       eq.topology.inside, profile, jphi);
            // psi_plasma must be the exact flux of the returned jphi, not the relaxed iterate.
            psi_pl = plasma_flux(jphi);
            eq.psi = total_psi();
            eq.psi_plasma = psi_pl;
            eq.jphi = jphi;
            eq.converged = true;
            return eq;
        }
    }

    eq.psi_plasma = psi_pl;
    eq.jphi = jphi;
    eq.converged = false;
    return eq;
}

} // namespace shafranov
