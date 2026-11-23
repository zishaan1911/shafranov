#pragma once

#include "core/equilibrium.h"
#include "core/machine.h"

#include <Eigen/Dense>

#include <vector>

namespace shafranov {

// Self-inductance (H) of a one-turn ring of major radius R with a rectangular dR x dZ
// cross-section: mu0 R (ln(8 R / g) - 2), with g = 0.2235 (dR + dZ) the rectangle's geometric mean
// distance from itself.
double self_inductance(double R, double dR, double dZ);

// Lumped circuit model of every conductor. Conductor order: the circuits of MachineGrid::circuits,
// then the vessel filaments. Each vessel filament is its own passive circuit.
struct CircuitModel {
    int n_circuits = 0;
    int n_vessel = 0;
    // Mutual inductances (H). Between distinct filaments: 2 pi greens_psi, weighted by turns; a
    // filament with itself: turns^2 * self_inductance of its own rectangle (a coil's sub-filament
    // rectangles are dR/nr x dZ/nz, as in circuit_windings).
    Eigen::MatrixXd M;
    // Resistance (ohm): a circuit's is the sum of its coils' resistances; a vessel filament's own.
    Eigen::VectorXd R;
};

CircuitModel make_circuit_model(const Machine& m, const MachineGrid& mg);

// Decay times (s) of the vessel's free current modes (L/R eigenvalues of the vessel block alone),
// longest first, with each mode's up-down parity: +1 for even (symmetric in Z), -1 for odd.
struct VesselMode {
    double tau;
    int parity;
};
std::vector<VesselMode> vessel_modes(const CircuitModel& c, const Machine& m);

} // namespace shafranov
