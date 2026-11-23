#include "core/circuits.h"

#include "core/coils.h"
#include "core/constants.h"
#include "core/greens.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace shafranov {

namespace {

struct Piece {
    double R, Z, turns, dR, dZ;
};

} // namespace

double self_inductance(double R, double dR, double dZ) {
    return mu0 * R * (std::log(8.0 * R / (0.2235 * (dR + dZ))) - 2.0);
}

CircuitModel make_circuit_model(const Machine& m, const MachineGrid& mg) {
    const std::size_t nc = mg.circuits.size();
    const std::size_t nv = m.vessel.size();
    const std::size_t n = nc + nv;

    std::vector<std::vector<Piece>> pieces(n);
    Eigen::VectorXd res = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(n));
    for (const Coil& c : m.coils) {
        const auto k = static_cast<std::size_t>(
            std::find(mg.circuits.begin(), mg.circuits.end(), c.circuit) - mg.circuits.begin());
        if (k >= nc) {
            throw std::runtime_error("coil " + c.name + " is in no circuit of the machine grid");
        }
        const int nr = std::max(1, static_cast<int>(std::ceil(c.dR / 0.02)));
        const int nz = std::max(1, static_cast<int>(std::ceil(c.dZ / 0.02)));
        const double turns = c.turns / (nr * nz);
        for (int a = 0; a < nr; ++a) {
            for (int b = 0; b < nz; ++b) {
                pieces[k].push_back({c.R + c.dR * ((a + 0.5) / nr - 0.5),
                                     c.Z + c.dZ * ((b + 0.5) / nz - 0.5), turns, c.dR / nr,
                                     c.dZ / nz});
            }
        }
        res(static_cast<Eigen::Index>(k)) += c.resistance;
    }
    for (std::size_t i = 0; i < nv; ++i) {
        const Filament& f = m.vessel[i];
        pieces[nc + i].push_back({f.R, f.Z, 1.0, f.dR, f.dZ});
        res(static_cast<Eigen::Index>(nc + i)) = f.resistance;
    }

    CircuitModel c;
    c.n_circuits = static_cast<int>(nc);
    c.n_vessel = static_cast<int>(nv);
    c.M = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
    c.R = res;

    // Workers take whole rows (upper triangle, mirrored); rows are uneven, so hand them out
    // one at a time.
    std::atomic<std::size_t> next{0};
    auto work = [&]() {
        for (std::size_t a = next++; a < n; a = next++) {
            for (std::size_t b = a; b < n; ++b) {
                double sum = 0.0;
                for (std::size_t i = 0; i < pieces[a].size(); ++i) {
                    const Piece& f = pieces[a][i];
                    for (std::size_t j = 0; j < pieces[b].size(); ++j) {
                        const Piece& h = pieces[b][j];
                        const double l = (a == b && i == j)
                                             ? self_inductance(f.R, f.dR, f.dZ)
                                             : mutual_inductance(f.R, f.Z, h.R, h.Z);
                        sum += f.turns * h.turns * l;
                    }
                }
                c.M(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = sum;
                c.M(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) = sum;
            }
        }
    };
    const unsigned workers = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> threads;
    for (unsigned w = 1; w < workers; ++w) {
        threads.emplace_back(work);
    }
    work();
    for (std::thread& t : threads) {
        t.join();
    }
    return c;
}

std::vector<VesselMode> vessel_modes(const CircuitModel& c, const Machine& m) {
    const Eigen::Index nv = c.n_vessel;
    if (static_cast<std::size_t>(nv) != m.vessel.size()) {
        throw std::runtime_error("circuit model and machine disagree on the vessel");
    }
    const Eigen::MatrixXd Lv = c.M.bottomRightCorner(nv, nv);
    const Eigen::VectorXd Rv = c.R.tail(nv);
    const Eigen::VectorXd d = Rv.cwiseSqrt().cwiseInverse();
    const Eigen::MatrixXd S = d.asDiagonal() * Lv * d.asDiagonal();
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("vessel eigenproblem did not converge");
    }

    std::vector<Eigen::Index> mirror(static_cast<std::size_t>(nv));
    for (Eigen::Index i = 0; i < nv; ++i) {
        const Filament& f = m.vessel[static_cast<std::size_t>(i)];
        Eigen::Index found = -1;
        for (Eigen::Index j = 0; j < nv; ++j) {
            const Filament& g = m.vessel[static_cast<std::size_t>(j)];
            if (std::abs(g.R - f.R) < 1e-6 && std::abs(g.Z + f.Z) < 1e-6) {
                found = j;
                break;
            }
        }
        if (found < 0) {
            throw std::runtime_error("vessel filament has no mirror image");
        }
        mirror[static_cast<std::size_t>(i)] = found;
    }

    std::vector<VesselMode> modes;
    for (Eigen::Index k = nv - 1; k >= 0; --k) {
        const Eigen::VectorXd x = d.cwiseProduct(solver.eigenvectors().col(k));
        double overlap = 0.0;
        for (Eigen::Index i = 0; i < nv; ++i) {
            overlap += x(i) * x(mirror[static_cast<std::size_t>(i)]);
        }
        modes.push_back({solver.eigenvalues()(k), overlap > 0.0 ? 1 : -1});
    }
    std::stable_sort(modes.begin(), modes.end(),
                     [](const VesselMode& a, const VesselMode& b) { return a.tau > b.tau; });
    return modes;
}

} // namespace shafranov
