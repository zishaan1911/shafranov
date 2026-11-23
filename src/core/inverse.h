#pragma once

#include "core/equilibrium.h"

#include <string>
#include <vector>

namespace shafranov {

// Shape targets for choosing coil currents, as in FreeGS's constraints (and TCV's FBT).
struct ShapeTargets {
    // Each pair of points should lie on the same flux surface: psi(R1, Z1) == psi(R2, Z2).
    struct Isoflux {
        double R1, Z1, R2, Z2;
    };
    std::vector<Isoflux> isoflux;
    std::vector<Point> xpoints;              // B_R = B_Z = 0 here
    std::vector<std::string> fixed;          // circuits left at their current value
    double gamma = 1e-12;                    // Tikhonov weight on the adjusted currents (A^-2)
};

// Feedback for solve_equilibrium that, every iteration, sets the free circuits' currents to the
// least-squares solution of the targets given the present plasma current (FreeGS's
// constrain step). Throws std::invalid_argument for unknown circuit names.
Feedback shape_feedback(const MachineGrid& mg, const ShapeTargets& targets);

} // namespace shafranov
