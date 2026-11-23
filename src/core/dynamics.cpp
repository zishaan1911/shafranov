#include "core/dynamics.h"

#include "core/constants.h"

#include <cstddef>
#include <string>
#include <utility>

namespace shafranov {

Eigen::VectorXd plasma_linkage(const MachineGrid& mg, const std::vector<double>& jphi) {
    const std::size_t n_circuits = mg.circuit_psi.size();
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(
        n_circuits + mg.vessel_psi.size()));
    const double dA = mg.grid.dR() * mg.grid.dZ();
    for (std::size_t k = 0; k < jphi.size(); ++k) {
        if (jphi[k] == 0.0) {
            continue;
        }
        const double w = 2.0 * pi * jphi[k] * dA;
        for (std::size_t c = 0; c < n_circuits; ++c) {
            phi[static_cast<Eigen::Index>(c)] += w * mg.circuit_psi[c][k];
        }
        for (std::size_t v = 0; v < mg.vessel_psi.size(); ++v) {
            phi[static_cast<Eigen::Index>(n_circuits + v)] += w * mg.vessel_psi[v][k];
        }
    }
    return phi;
}

Simulation::Simulation(const MachineGrid& mg, const CircuitModel& model, const Profile& profile,
                       const Equilibrium& start, const Currents& start_currents, double dt)
    : mg_(mg), model_(model), profile_(profile), eq_(start), currents_(start_currents), dt_(dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument("Simulation: dt must be positive");
    }
    if (currents_.vessel.empty()) {
        currents_.vessel.assign(static_cast<std::size_t>(model.n_vessel), 0.0);
    }
    if (currents_.circuits.size() != static_cast<std::size_t>(model.n_circuits) ||
        currents_.vessel.size() != static_cast<std::size_t>(model.n_vessel) ||
        mg.circuit_psi.size() != static_cast<std::size_t>(model.n_circuits) ||
        mg.vessel_psi.size() != static_cast<std::size_t>(model.n_vessel) ||
        eq_.jphi.size() != static_cast<std::size_t>(mg.grid.size())) {
        throw std::invalid_argument("Simulation: size mismatch between grid, model and currents");
    }
    phi_ = plasma_linkage(mg_, eq_.jphi);
    lhs(0); // fail early if M + dt R is not positive definite
}

const Eigen::LLT<Eigen::MatrixXd>& Simulation::lhs(int depth) {
    while (static_cast<int>(lhs_.size()) <= depth) {
        const double h = dt_ / static_cast<double>(1 << lhs_.size());
        Eigen::MatrixXd a = model_.M;
        a.diagonal() += h * model_.R;
        lhs_.emplace_back(a);
        if (lhs_.back().info() != Eigen::Success) {
            throw std::runtime_error("Simulation: M + dt R is not positive definite");
        }
    }
    return lhs_[static_cast<std::size_t>(depth)];
}

std::vector<double> Simulation::holding_voltages() const {
    std::vector<double> v(static_cast<std::size_t>(model_.n_circuits));
    for (std::size_t c = 0; c < v.size(); ++c) {
        v[c] = model_.R[static_cast<Eigen::Index>(c)] * currents_.circuits[c];
    }
    return v;
}

void Simulation::step(const std::vector<double>& voltages) {
    if (voltages.size() != static_cast<std::size_t>(model_.n_circuits)) {
        throw std::invalid_argument("Simulation::step needs one voltage per circuit");
    }
    substep(voltages, 0);
}

void Simulation::substep(const std::vector<double>& voltages, int depth) {
    std::string why;
    if (advance(voltages, depth, why)) {
        return;
    }
    if (depth == 2) {
        throw PlasmaLost(why);
    }
    substep(voltages, depth + 1);
    substep(voltages, depth + 1);
}

bool Simulation::advance(const std::vector<double>& voltages, int depth, std::string& why) {
    const double h = dt_ / static_cast<double>(1 << depth);
    const Eigen::LLT<Eigen::MatrixXd>& factor = lhs(depth);
    const std::size_t nc = static_cast<std::size_t>(model_.n_circuits);
    const std::size_t nv = static_cast<std::size_t>(model_.n_vessel);
    Eigen::VectorXd i_old(static_cast<Eigen::Index>(nc + nv));
    for (std::size_t c = 0; c < nc; ++c) {
        i_old[static_cast<Eigen::Index>(c)] = currents_.circuits[c];
    }
    for (std::size_t v = 0; v < nv; ++v) {
        i_old[static_cast<Eigen::Index>(nc + v)] = currents_.vessel[v];
    }
    const Eigen::VectorXd phi_old = phi_;

    Eigen::VectorXd drive = Eigen::VectorXd::Zero(i_old.size());
    for (std::size_t c = 0; c < nc; ++c) {
        drive[static_cast<Eigen::Index>(c)] = h * voltages[c];
    }
    const Eigen::VectorXd base = model_.M * i_old + drive + phi_old;

    const Feedback feedback = [&](const Equilibrium& eq, Currents& currents) {
        const Eigen::VectorXd phi = plasma_linkage(mg_, eq.jphi);
        const Eigen::VectorXd i_new = factor.solve(base - phi);
        for (std::size_t c = 0; c < nc; ++c) {
            currents.circuits[c] = i_new[static_cast<Eigen::Index>(c)];
        }
        for (std::size_t v = 0; v < nv; ++v) {
            currents.vessel[v] = i_new[static_cast<Eigen::Index>(nc + v)];
        }
    };

    Currents currents = currents_;
    const CriticalPoint& axis = eq_.topology.axis;
    InitialPlasma initial{axis.R, axis.Z, 0.1};
    initial.jphi = eq_.jphi;
    SolveOptions options;
    options.max_iterations = 300;
    // The plasma boundary moving across grid nodes leaves a ~1e-6 noise floor in the Picard
    // residual; a step's real change is ~1e-4 and above, so 1e-5 resolves it without chasing noise.
    options.tolerance = 1e-5;
    // A warm-started step changes the state only slightly: plain Picard at full relaxation
    // converges in ~15 iterations here, where Anderson extrapolation churns (~60).
    options.relax = 1.0;
    options.anderson = 0;

    const std::string when = " at t = " + std::to_string(time_ + h) + " s";
    Equilibrium next;
    try {
        next = solve_equilibrium(mg_, currents, profile_, initial, options, feedback);
    } catch (const std::runtime_error& e) {
        why = std::string("plasma lost") + when + ": " + e.what();
        return false;
    }
    if (!next.converged) {
        why = "the solver can no longer follow the plasma (axis at Z = " +
              std::to_string(eq_.topology.axis.Z) + " m)" + when;
        return false;
    }

    // Close the flux books exactly with the converged plasma: the iteration's last circuit solve
    // used the previous iterate's current, which differs by up to the tolerance, and that error
    // would otherwise accumulate from step to step.
    phi_ = plasma_linkage(mg_, next.jphi);
    const Eigen::VectorXd i_new = factor.solve(base - phi_);
    for (std::size_t c = 0; c < nc; ++c) {
        currents.circuits[c] = i_new[static_cast<Eigen::Index>(c)];
    }
    for (std::size_t v = 0; v < nv; ++v) {
        currents.vessel[v] = i_new[static_cast<Eigen::Index>(nc + v)];
    }
    eq_ = std::move(next);
    currents_ = std::move(currents);
    time_ += h;
    return true;
}

} // namespace shafranov
