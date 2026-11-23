#include "core/dynamics.h"

#include "core/eqdsk.h"
#include "sim/worker.h"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
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

// Slope (1/s) of log|Z - Z0| against t over the samples between 1 mm and 1 cm of displacement.
double growth_of(const std::vector<double>& t, const std::vector<double>& dz) {
    double n = 0.0, st = 0.0, sy = 0.0, stt = 0.0, sty = 0.0;
    for (std::size_t k = 0; k < t.size(); ++k) {
        const double d = std::abs(dz[k]);
        if (d > 1e-3 && d < 1e-2) {
            const double y = std::log(d);
            n += 1.0;
            st += t[k];
            sy += y;
            stt += t[k] * t[k];
            sty += t[k] * y;
        }
    }
    const double den = n * stt - st * st;
    if (n < 2.0 || den == 0.0) {
        return 0.0;
    }
    return (n * sty - st * sy) / den;
}

// Vertical growth rate after a +20 V kick on circuit G with every other voltage held.
double growth(double kappa) {
    Request r;
    r.kappa = kappa;
    const Start s = solve_start(r);
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    const double z0 = s.eq.topology.axis.Z;
    const std::vector<double> hold = sim.holding_voltages();
    std::vector<double> kick = hold;
    kick[circuit_index("G")] += 20.0;

    std::vector<double> t, dz;
    for (int k = 0; k < 400; ++k) {
        try {
            sim.step(k == 0 ? kick : hold);
        } catch (const PlasmaLost&) {
            break;
        }
        const double d = sim.equilibrium().topology.axis.Z - z0;
        t.push_back(sim.time());
        dz.push_back(d);
        if (std::abs(d) > 0.02) {
            break;
        }
    }
    const double rate = growth_of(t, dz);
    MESSAGE("kappa = " << kappa << ": vertical growth rate " << rate << " 1/s");
    return rate;
}

} // namespace

TEST_CASE("a circular limited plasma held by the holding voltages stays put") {
    Request r;
    r.kappa = 1.0;
    r.diverted = false;
    const Start s = solve_start(r);
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    const double z0 = s.eq.topology.axis.Z;
    const double r0 = s.eq.topology.axis.R;
    const std::vector<double> hold = sim.holding_voltages();
    for (int k = 0; k < 20; ++k) {
        sim.step(hold);
    }
    const CriticalPoint& axis = sim.equilibrium().topology.axis;
    CHECK(std::abs(axis.Z - z0) < 1e-4);
    CHECK(std::abs(axis.R - r0) < 1e-4);
    for (const double i : sim.currents().vessel) {
        CHECK(std::abs(i) < 1.0);
    }
}

TEST_CASE("an elongated plasma is vertically unstable") {
    const double rate = growth(1.4);
    CHECK(rate > 50.0);
    CHECK(rate < 5000.0);
}

TEST_CASE("elongation destabilises the vertical position") {
    // 1.5, not higher: by kappa 1.6 the growth time nears dt = 0.1 ms and backward Euler is too coarse.
    CHECK(growth(1.5) > 1.5 * growth(1.3));
}

TEST_CASE("a step conserves flux linkage") {
    const Start s = solve_start(Request{});
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    const CircuitModel& model = tcv_model();
    const std::size_t nc = static_cast<std::size_t>(model.n_circuits);
    const std::size_t nv = static_cast<std::size_t>(model.n_vessel);

    const auto gather = [&](const Currents& c) {
        Eigen::VectorXd i(static_cast<Eigen::Index>(nc + nv));
        for (std::size_t k = 0; k < nc; ++k) {
            i[static_cast<Eigen::Index>(k)] = c.circuits[k];
        }
        for (std::size_t k = 0; k < nv; ++k) {
            i[static_cast<Eigen::Index>(nc + k)] = c.vessel[k];
        }
        return i;
    };

    std::vector<double> v = sim.holding_voltages();
    v[circuit_index("G")] += 20.0;
    const Eigen::VectorXd i_old = gather(sim.currents());
    const Eigen::VectorXd phi_old = plasma_linkage(tcv_grid(), sim.equilibrium().jphi);
    sim.step(v);
    const Eigen::VectorXd i_new = gather(sim.currents());
    const Eigen::VectorXd phi_new = plasma_linkage(tcv_grid(), sim.equilibrium().jphi);

    Eigen::VectorXd volts = Eigen::VectorXd::Zero(i_old.size());
    for (std::size_t k = 0; k < nc; ++k) {
        volts[static_cast<Eigen::Index>(k)] = v[k];
    }
    const Eigen::VectorXd residual =
        model.M * (i_new - i_old) + (phi_new - phi_old) + dt * model.R.cwiseProduct(i_new) -
        dt * volts;
    CHECK(residual.cwiseAbs().maxCoeff() < 1e-6 * (model.M * i_old).cwiseAbs().maxCoeff());
}

TEST_CASE("a vertical displacement event ends with a lost plasma") {
    Request r;
    r.kappa = 1.6;
    const Start s = solve_start(r);
    Simulation sim(tcv_grid(), tcv_model(), s.profile, s.eq, s.currents, dt);
    const std::vector<double> hold = sim.holding_voltages();
    std::vector<double> kick = hold;
    kick[circuit_index("G")] += 20.0;

    bool lost = false;
    for (int k = 0; k < 3000 && !lost; ++k) {
        try {
            sim.step(k == 0 ? kick : hold);
        } catch (const PlasmaLost&) {
            lost = true;
        }
    }
    CHECK(lost);
}
