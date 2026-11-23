#include "sim/discharge.h"

#include "core/eqdsk.h"

#include <doctest/doctest.h>

#include <chrono>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace shafranov;

namespace {

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

std::vector<Waypoint> short_scenario() {
    Request a;
    a.Ip = 150e3;
    a.kappa = 1.3;
    Request b;
    b.Ip = 160e3;
    b.kappa = 1.32;
    return {{0.0, a}, {0.002, b}};
}

// Polls until the view's status is one of the two, or the timeout passes.
std::shared_ptr<const DischargeView> wait_for(const Discharge& d, DischargeStatus a,
                                              DischargeStatus b, int seconds, int poll_ms) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(seconds);
    while (std::chrono::steady_clock::now() < deadline) {
        const auto v = d.view();
        if (v->status == a || v->status == b) {
            return v;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(poll_ms));
    }
    return d.view();
}

} // namespace

TEST_CASE("make_plan solves the demo scenario and interpolates its references") {
    const Plan plan = make_plan(tcv_grid(), demo_scenario());
    REQUIRE(plan.currents.size() == 3);
    CHECK(plan.start.converged);
    CHECK(reference_ip(plan, 0.01) == doctest::Approx(160e3).epsilon(1e-12));
    CHECK(reference_currents(plan, -1.0) == plan.currents.front());
    CHECK(reference_currents(plan, 1.0) == plan.currents.back());
}

TEST_CASE("make_plan rejects bad waypoint lists") {
    CHECK_THROWS_AS(make_plan(tcv_grid(), {}), std::invalid_argument);

    std::vector<Waypoint> same = short_scenario();
    same[1].t = same[0].t;
    CHECK_THROWS_AS(make_plan(tcv_grid(), same), std::invalid_argument);

    std::vector<Waypoint> backwards = short_scenario();
    backwards[1].t = -0.001;
    CHECK_THROWS_AS(make_plan(tcv_grid(), backwards), std::invalid_argument);
}

TEST_CASE("a short discharge runs to the end and stays vertically stable") {
    Discharge d(tcv_grid(), tcv_model(), short_scenario(), 0.003);
    const double z0 = d.view()->eq.topology.axis.Z;
    d.run();
    const auto v = wait_for(d, DischargeStatus::Finished, DischargeStatus::Lost, 120, 20);
    REQUIRE(v->status == DischargeStatus::Finished);
    CHECK(static_cast<int>(v->history.size()) >= 29);
    CHECK(static_cast<int>(v->history.size()) <= 31);
    for (const Sample& s : v->history) {
        CHECK(std::abs(s.Z_axis - z0) < 2e-3);
    }
    CHECK(v->history.back().Ip == doctest::Approx(160e3).epsilon(1e-9));
}

TEST_CASE("reset returns a finished discharge to the start") {
    Discharge d(tcv_grid(), tcv_model(), short_scenario(), 0.003);
    d.run();
    const auto done = wait_for(d, DischargeStatus::Finished, DischargeStatus::Lost, 120, 20);
    REQUIRE(done->status == DischargeStatus::Finished);

    d.reset();
    const auto v = wait_for(d, DischargeStatus::Paused, DischargeStatus::Paused, 5, 20);
    CHECK(v->status == DischargeStatus::Paused);
    CHECK(v->t == 0.0);
    CHECK(v->history.empty());
}

TEST_CASE("the vertical target follows the planned axis height") {
    const Plan p = make_plan(tcv_grid(), demo_scenario());
    REQUIRE(p.axis_z.size() == p.waypoints.size());
    CHECK(reference_z(p, -1.0) == p.axis_z.front());
    CHECK(reference_z(p, 0.01) == doctest::Approx(0.5 * (p.axis_z[0] + p.axis_z[1])));
    CHECK(p.axis_z[0] == p.start.topology.axis.Z);
}
