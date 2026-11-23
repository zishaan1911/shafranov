#pragma once

#include "core/coils.h"
#include "core/grid.h"
#include "core/machine.h"
#include "core/profiles.h"
#include "core/topology.h"

#include <functional>
#include <string>
#include <vector>

namespace shafranov {

// Flux maps that depend only on the machine and the grid. Expensive to build; reuse it.
struct MachineGrid {
    Grid grid;
    std::vector<std::string> circuits;              // circuit_names(machine)
    std::vector<std::vector<Winding>> windings;     // circuit_windings(machine)
    std::vector<std::vector<double>> circuit_psi;   // psi on the grid per 1 A in each circuit
    std::vector<std::vector<double>> vessel_psi;    // psi on the grid per 1 A in each vessel filament
    std::vector<int> boundary;                      // grid indices of the boundary nodes
    std::vector<double> boundary_greens;            // [b * grid.size() + k]: psi at boundary node b
                                                    // from 1 A at grid node k
    std::vector<Point> limiter;
    std::vector<Filament> vessel;
};

MachineGrid make_machine_grid(const Machine& m, const Grid& g);

// Conductor currents (A): one per circuit (MachineGrid::circuits order) and one per vessel
// filament (empty means all zero).
struct Currents {
    std::vector<double> circuits;
    std::vector<double> vessel;
};

struct Equilibrium {
    std::vector<double> psi;        // total poloidal flux per radian on the grid
    std::vector<double> psi_plasma; // the plasma's own contribution
    std::vector<double> jphi;       // toroidal current density (A/m^2)
    Topology topology;
    double lambda = 0.0;            // profile scale (see plasma_current)
    int iterations = 0;
    bool converged = false;
};

// First guess for the plasma current: a parabolic disc of radius a centred at (R, Z), or, if jphi
// is not empty, that current density (e.g. a previous solution, for a warm start) rescaled to Ip.
struct InitialPlasma {
    double R, Z, a;
    std::vector<double> jphi = {};
};

struct SolveOptions {
    int max_iterations = 500; // weakly restoring fields converge slowly (contraction ~0.9)
    double tolerance = 1e-6; // on max |change in psi| / |psi_axis - psi_bndry|
    double relax = 0.5;      // weight of the new plasma flux in each Picard update
};

// Called once per iteration with the latest equilibrium; may change the circuit currents (shape
// and position control). Leave empty to hold the currents fixed.
using Feedback = std::function<void(const Equilibrium&, Currents&)>;

// Free-boundary Picard iteration. Throws std::invalid_argument on size mismatches and
// std::runtime_error if the plasma is lost (no axis inside the limiter).
Equilibrium solve_equilibrium(const MachineGrid& mg, Currents& currents, const Profile& profile,
                              const InitialPlasma& initial, const SolveOptions& options = {},
                              const Feedback& feedback = {});

} // namespace shafranov
