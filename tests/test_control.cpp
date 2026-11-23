#include "core/control.h"

#include "core/dynamics.h"
#include "core/eqdsk.h"
#include "sim/worker.h"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace shafranov;

namespace {

constexpr double dt = 1e-4;

// Slow to build, so shared by every test case in this file.
const Machine& tcv_machine() {
    static const Machine m = load_machine(SHAFRANOV_DATA_DIR "/tcv");
    return m;
}

const MachineGrid& tcv_grid() {
    static const MachineGrid mg =
        make_machine_grid(tcv_machine(), Grid(65, 129, 0.55, 1.20, -0.80, 0.80));
    return mg;
}

const CircuitModel& tcv_model() {
    static const CircuitModel model = make_circuit_model(tcv_machine(), tcv_grid());
    return model;
}

struct Start {
    Profile profile;
    Equilibrium eq;
    Currents currents;
};

// The inverse solve for the request's shape, from a cold start, vessel currents zero.
Start solve_start(const Request& r) {
    const MachineGrid& mg = tcv_grid();
    Start s;
    s.profile = Profile{r.Ip, r.beta0, 1.0, 2.0, r.R0};
    s.currents.circuits.assign(mg.circuits.size(), 0.0);
    const InitialPlasma initial{r.R0, 0.0, 0.7 * r.a};
    s.eq = solve_equilibrium(mg, s.currents, s.profile, initial, SolveOptions{},
                             shape_feedback(mg, targets_for(r)));
    s.currents.vessel.assign(static_cast<std::size_t>(tcv_model().n_vessel), 0.0);
    return s;
}

std::size_t circuit_index(const std::string& name) {
    const auto& names = tcv_grid().circuits;
    return static_cast<std::size_t>(std::find(names.begin(), names.end(), name) - names.begin());
}

// Synthetic circuit currents (A), small enough that no voltage reaches its supply limit.
Currents synthetic_currents() {
    Currents c;
    c.circuits.resize(tcv_grid().circuits.size());
    for (std::size_t k = 0; k < c.circuits.size(); ++k) {
        c.circuits[k] = 10.0 * static_cast<double>(k + 1);
    }
    return c;
}

// Runs the controller on a kicked kappa 1.5 plasma for up to `steps` steps; the first step
// carries `kick` extra volts on G. Returns the axis height error at the end.
double run_kicked(Controller& ctrl, int steps, double kick, std::vector<double> reference,
                  const Start& s, Simulation& sim) {
    const double z0 = s.eq.topology.axis.Z;
    ctrl.vertical.z_ref = z0;
    double z_prev = z0;
    for (int k = 0; k < steps; ++k) {
        const double z = sim.equilibrium().topology.axis.Z;
        std::vector<double> v = ctrl.voltages(sim.currents(), reference, z, z_prev, dt);
        if (k == 0) {
            v[circuit_index("G")] += kick;
        }
        z_prev = z;
        sim.step(v);
    }
    return sim.equilibrium().topology.axis.Z - z0;
}

} // namespace

TEST_CASE("supply limits follow the circuit names") {
    const std::vector<double> limits = supply_limits({"OH1", "OH2", "E1", "F8", "G", "T"});
    const std::vector<double> expected = {1400.0, 1400.0, 651.0, 1250.0, 566.0, 0.0};
    REQUIRE(limits.size() == expected.size());
    for (std::size_t k = 0; k < expected.size(); ++k) {
        CHECK(limits[k] == doctest::Approx(expected[k]));
    }
}

TEST_CASE("with the plant at its reference the voltages are R I") {
    const CircuitModel& model = tcv_model();
    const std::vector<std::string>& names = tcv_grid().circuits;
    const Controller ctrl(model, names);
    const Currents now = synthetic_currents();
    const double z_ref = ctrl.vertical.z_ref;

    const std::vector<double> v = ctrl.voltages(now, now.circuits, z_ref, z_ref, dt);
    const std::vector<double> limits = supply_limits(names);
    REQUIRE(v.size() == names.size());
    for (std::size_t k = 0; k < v.size(); ++k) {
        if (names[k] == "T") {
            CHECK(v[k] == 0.0);
        } else if (limits[k] > 0.0) {
            const double expected = model.R[static_cast<Eigen::Index>(k)] * now.circuits[k];
            CHECK(v[k] == doctest::Approx(expected).epsilon(1e-12));
        }
    }

    const std::size_t e1 = circuit_index("E1");
    std::vector<double> reference = now.circuits;
    reference[e1] += 10.0;
    const std::vector<double> v2 = ctrl.voltages(now, reference, z_ref, z_ref, dt);
    const auto i = static_cast<Eigen::Index>(e1);
    const double change = 10.0 * (model.R[i] + model.M(i, i) / ctrl.current.tau);
    CHECK(v2[e1] - v[e1] == doctest::Approx(change).epsilon(1e-12));
}

TEST_CASE("voltages are clipped to the supply limits") {
    const Controller ctrl(tcv_model(), tcv_grid().circuits);
    const Currents now = synthetic_currents();
    const std::size_t e1 = circuit_index("E1");
    const double z_ref = ctrl.vertical.z_ref;

    std::vector<double> reference = now.circuits;
    reference[e1] += 1e6;
    CHECK(ctrl.voltages(now, reference, z_ref, z_ref, dt)[e1] == 651.0);
    reference[e1] = now.circuits[e1] - 1e6;
    CHECK(ctrl.voltages(now, reference, z_ref, z_ref, dt)[e1] == -651.0);
}

TEST_CASE("vertical control holds a kicked elongated plasma") {
    Request r;
    r.kappa = 1.5;
    const Start s = solve_start(r);
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    Controller ctrl(tcv_model(), tcv_grid().circuits);

    double dz = 0.0;
    CHECK_NOTHROW(dz = run_kicked(ctrl, 100, 20.0, s.currents.circuits, s, sim));
    CHECK(std::abs(dz) < 1e-3);
}

TEST_CASE("without vertical control the same kick loses the plasma") {
    Request r;
    r.kappa = 1.5;
    const Start s = solve_start(r);
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    Controller ctrl(tcv_model(), tcv_grid().circuits);
    ctrl.vertical_on = false;

    CHECK_THROWS_AS(run_kicked(ctrl, 300, 20.0, s.currents.circuits, s, sim), PlasmaLost);
}

TEST_CASE("coil currents track their reference") {
    Request r;
    r.kappa = 1.5;
    const Start s = solve_start(r);
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    Controller ctrl(tcv_model(), tcv_grid().circuits);
    ctrl.current.tau = 2e-3;

    const std::size_t e3 = circuit_index("E3");
    std::vector<double> reference = s.currents.circuits;
    reference[e3] += 200.0;
    run_kicked(ctrl, 60, 0.0, reference, s, sim);

    const double error = std::abs(sim.currents().circuits[e3] - reference[e3]);
    MESSAGE("E3 error after 60 steps: " << error << " A");
    // 60 steps are 3 tau: an ideal first-order lag leaves e^-3 = 5% (10 A); vessel coupling slows it.
    CHECK(error < 20.0);
}
