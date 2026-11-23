#include "core/inverse.h"

#include "core/greens.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace shafranov {

namespace {

enum class Quantity { Psi, BR, BZ };

struct Term {
    double R, Z, sign;
};

// One least-squares row: sum over terms of sign * quantity(term) == 0.
struct Row {
    Quantity quantity;
    std::vector<Term> terms;
    std::vector<double> circuit; // per amp in each circuit
    std::vector<double> vessel;  // per amp in each vessel filament
};

double filament_value(Quantity q, double Rc, double Zc, double R, double Z) {
    if (q == Quantity::Psi) {
        return greens_psi(Rc, Zc, R, Z);
    }
    const BField b = greens_b(Rc, Zc, R, Z);
    return q == Quantity::BR ? b.R : b.Z;
}

Row make_row(const MachineGrid& mg, Quantity q, std::vector<Term> terms) {
    Row row;
    row.quantity = q;
    row.terms = std::move(terms);
    row.circuit.assign(mg.windings.size(), 0.0);
    for (const Term& t : row.terms) {
        for (std::size_t c = 0; c < mg.windings.size(); ++c) {
            for (const Winding& w : mg.windings[c]) {
                row.circuit[c] += t.sign * w.turns * filament_value(q, w.R, w.Z, t.R, t.Z);
            }
        }
    }
    row.vessel.assign(mg.vessel.size(), 0.0);
    for (const Term& t : row.terms) {
        for (std::size_t v = 0; v < mg.vessel.size(); ++v) {
            row.vessel[v] += t.sign * filament_value(q, mg.vessel[v].R, mg.vessel[v].Z, t.R, t.Z);
        }
    }
    return row;
}

// Solves the n x n system M x = rhs in place by Gaussian elimination with partial pivoting.
std::vector<double> solve_dense(std::vector<double> M, std::vector<double> rhs, std::size_t n) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < n; ++r) {
            if (std::abs(M[r * n + col]) > std::abs(M[pivot * n + col])) {
                pivot = r;
            }
        }
        if (M[pivot * n + col] == 0.0) {
            throw std::runtime_error("shape_feedback: singular normal equations");
        }
        if (pivot != col) {
            for (std::size_t k = 0; k < n; ++k) {
                std::swap(M[pivot * n + k], M[col * n + k]);
            }
            std::swap(rhs[pivot], rhs[col]);
        }
        for (std::size_t r = col + 1; r < n; ++r) {
            const double f = M[r * n + col] / M[col * n + col];
            for (std::size_t k = col; k < n; ++k) {
                M[r * n + k] -= f * M[col * n + k];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    std::vector<double> x(n);
    for (std::size_t i = n; i-- > 0;) {
        double s = rhs[i];
        for (std::size_t k = i + 1; k < n; ++k) {
            s -= M[i * n + k] * x[k];
        }
        x[i] = s / M[i * n + i];
    }
    return x;
}

} // namespace

Feedback shape_feedback(const MachineGrid& mg, const ShapeTargets& targets) {
    const std::size_t nc = mg.circuits.size();
    std::vector<char> is_fixed(nc, 0);
    for (const std::string& name : targets.fixed) {
        const auto it = std::find(mg.circuits.begin(), mg.circuits.end(), name);
        if (it == mg.circuits.end()) {
            throw std::invalid_argument("shape_feedback: unknown circuit " + name);
        }
        is_fixed[static_cast<std::size_t>(it - mg.circuits.begin())] = 1;
    }
    std::vector<std::size_t> free_idx;
    for (std::size_t c = 0; c < nc; ++c) {
        if (!is_fixed[c]) {
            free_idx.push_back(c);
        }
    }

    std::vector<Row> rows;
    for (const auto& p : targets.isoflux) {
        rows.push_back(make_row(mg, Quantity::Psi, {{p.R1, p.Z1, 1.0}, {p.R2, p.Z2, -1.0}}));
    }
    for (const Point& x : targets.xpoints) {
        rows.push_back(make_row(mg, Quantity::BR, {{x.R, x.Z, 1.0}}));
        rows.push_back(make_row(mg, Quantity::BZ, {{x.R, x.Z, 1.0}}));
    }
    const double gamma = targets.gamma;

    return [&mg, rows, free_idx, is_fixed, gamma](const Equilibrium& eq, Currents& currents) {
        const std::size_t n = free_idx.size();
        const std::size_t m = rows.size();
        if (n == 0 || m == 0) {
            return;
        }
        const Grid& g = mg.grid;
        const double cell = g.dR() * g.dZ();

        std::vector<double> A(m * n);
        std::vector<double> b(m);
        for (std::size_t r = 0; r < m; ++r) {
            const Row& row = rows[r];
            double rest = 0.0;
            for (const Term& t : row.terms) {
                double plasma = 0.0;
                for (int i = 0; i < g.nR; ++i) {
                    for (int j = 0; j < g.nZ; ++j) {
                        const double J = eq.jphi[static_cast<std::size_t>(g.idx(i, j))];
                        if (J != 0.0) {
                            plasma += filament_value(row.quantity, g.R(i), g.Z(j), t.R, t.Z) * J *
                                      cell;
                        }
                    }
                }
                rest += t.sign * plasma;
            }
            for (std::size_t c = 0; c < is_fixed.size(); ++c) {
                if (is_fixed[c]) {
                    rest += currents.circuits[c] * row.circuit[c];
                }
            }
            for (std::size_t v = 0; v < currents.vessel.size(); ++v) {
                rest += currents.vessel[v] * row.vessel[v];
            }
            for (std::size_t k = 0; k < n; ++k) {
                A[r * n + k] = row.circuit[free_idx[k]];
            }
            b[r] = -rest;
        }

        std::vector<double> N(n * n, 0.0);
        std::vector<double> rhs(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t r = 0; r < m; ++r) {
                    N[i * n + j] += A[r * n + i] * A[r * n + j];
                }
            }
            N[i * n + i] += gamma;
            for (std::size_t r = 0; r < m; ++r) {
                rhs[i] += A[r * n + i] * b[r];
            }
        }
        const std::vector<double> x = solve_dense(std::move(N), std::move(rhs), n);
        for (std::size_t k = 0; k < n; ++k) {
            currents.circuits[free_idx[k]] = x[k];
        }
    };
}

} // namespace shafranov
