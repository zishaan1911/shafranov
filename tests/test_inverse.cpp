#include "core/inverse.h"

#include "core/constants.h"

#include <doctest/doctest.h>

#include <cmath>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

using namespace shafranov;

namespace {

// Slow to build, so shared by every test case in this file.
const MachineGrid& tcv_grid() {
    static const MachineGrid mg = make_machine_grid(load_machine(SHAFRANOV_DATA_DIR "/tcv"),
                                                    Grid(65, 129, 0.55, 1.20, -0.80, 0.80));
    return mg;
}

const Point xtarget{0.814, -0.339};

std::vector<Point> boundary_targets() {
    std::vector<Point> pts;
    for (int k = 0; k < 12; ++k) {
        if (k == 9) {
            continue; // the bottom, replaced by the X-point
        }
        const double t = k * 2.0 * pi / 12.0;
        pts.push_back({0.88 + 0.22 * std::cos(t + 0.3 * std::sin(t)), 1.4 * 0.22 * std::sin(t)});
    }
    return pts;
}

ShapeTargets tcv_targets() {
    ShapeTargets targets;
    for (const Point& p : boundary_targets()) {
        targets.isoflux.push_back({xtarget.R, xtarget.Z, p.R, p.Z});
    }
    targets.xpoints.push_back(xtarget);
    targets.fixed = {"OH1", "OH2", "G", "T"};
    targets.gamma = 1e-12;
    return targets;
}

struct Solved {
    Currents currents;
    Equilibrium eq;
};

// The lower single-null case of tools/freegs_reference.py, solved once for every test case.
const Solved& tcv_lsn() {
    static const Solved s = [] {
        const MachineGrid& mg = tcv_grid();
        Solved out;
        out.currents.circuits.assign(mg.circuits.size(), 0.0);
        const Profile profile{2e5, 0.3, 1.0, 2.0, 0.88};
        const InitialPlasma initial{0.88, 0.0, 0.15};
        out.eq = solve_equilibrium(mg, out.currents, profile, initial, SolveOptions{},
                                   shape_feedback(mg, tcv_targets()));
        return out;
    }();
    return s;
}

} // namespace

TEST_CASE("shape_feedback finds currents for a lower single-null TCV plasma") {
    const MachineGrid& mg = tcv_grid();
    const Currents& currents = tcv_lsn().currents;
    const Equilibrium& eq = tcv_lsn().eq;
    CHECK(eq.converged);

    bool near_target = false;
    for (const CriticalPoint& x : eq.topology.xpoints) {
        if (std::hypot(x.R - xtarget.R, x.Z - xtarget.Z) < 0.01) {
            near_target = true;
        }
    }
    CHECK(near_target);
    CHECK_FALSE(eq.topology.limited);
    CHECK(std::abs(eq.topology.axis.R - 0.88) < 0.05);
    CHECK(std::abs(eq.topology.axis.Z) < 0.05);

    const double psi_axis = eq.topology.axis.psi;
    for (const Point& p : boundary_targets()) {
        const double psi = interpolate(mg.grid, eq.psi, p.R, p.Z);
        const double psi_n = (psi - psi_axis) / (eq.topology.psi_bndry - psi_axis);
        // Regularisation trades shape for realistic currents: FreeGS misses these by up to 0.143.
        CHECK(std::abs(psi_n - 1.0) < 0.15);
    }

    const std::vector<std::string> fixed = {"OH1", "OH2", "G", "T"};
    for (std::size_t c = 0; c < mg.circuits.size(); ++c) {
        bool is_fixed = false;
        for (const std::string& name : fixed) {
            is_fixed = is_fixed || mg.circuits[c] == name;
        }
        if (is_fixed) {
            CHECK(currents.circuits[c] == 0.0);
        } else {
            CHECK(std::isfinite(currents.circuits[c]));
            CHECK(std::abs(currents.circuits[c]) < 7.7e3); // TCV E and F supply limit
        }
    }
}

TEST_CASE("shape_feedback rejects an unknown fixed circuit") {
    ShapeTargets targets;
    targets.fixed = {"NOPE"};
    CHECK_THROWS_AS(shape_feedback(tcv_grid(), targets), std::invalid_argument);
}

