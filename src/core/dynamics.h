#pragma once

#include "core/circuits.h"
#include "core/equilibrium.h"

#include <stdexcept>
#include <vector>

namespace shafranov {

// The plasma left the limiter or lost its magnetic axis (a vertical displacement event ends here).
struct PlasmaLost : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Time evolution of a massless plasma coupled to every conductor. Each step solves, by backward Euler,
//   (M + dt R) I_new = M I_old + dt V - (Phi_new - Phi_old)
// where Phi is the flux the plasma current links with each conductor (2 pi sum_k psi_c[k] J_k dA, by
// reciprocity from MachineGrid's flux maps), together with the plasma equilibrium for I_new. The
// circuit solve runs as solve_equilibrium's per-iteration Feedback, warm-started from the last step.
// The plasma current Ip is prescribed by the Profile (it does not evolve).
class Simulation {
public:
    // start: a converged equilibrium for start_currents (circuits, and vessel: empty = zeros).
    Simulation(const MachineGrid& mg, const CircuitModel& model, const Profile& profile,
               const Equilibrium& start, const Currents& start_currents, double dt);

    // Advances one step of dt with the given circuit voltages (size n_circuits). Throws PlasmaLost
    // when the plasma is lost (state() then still holds the last good step).
    void step(const std::vector<double>& voltages);
    // Voltages that hold the starting circuit currents in steady state: V = R I.
    std::vector<double> holding_voltages() const;

    double time() const { return time_; }
    const Equilibrium& equilibrium() const { return eq_; }
    const Currents& currents() const { return currents_; }
    double dt() const { return dt_; }
    void set_profile(const Profile& p) { profile_ = p; }

private:
    const MachineGrid& mg_;
    const CircuitModel& model_;
    Profile profile_;
    Equilibrium eq_;
    Currents currents_;
    Eigen::VectorXd phi_; // plasma flux linkage with each conductor at the current state
    Eigen::LLT<Eigen::MatrixXd> lhs_;
    double dt_;
    double time_ = 0.0;
};

// Plasma flux linkage (Wb) with each conductor, in CircuitModel order.
Eigen::VectorXd plasma_linkage(const MachineGrid& mg, const std::vector<double>& jphi);

} // namespace shafranov
