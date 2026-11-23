#include "sim/worker.h"

#include "core/eqdsk.h"

#include <doctest/doctest.h>

#include <chrono>
#include <cmath>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace shafranov;

namespace {

// Slow to build, so shared by every test case in this file.
const MachineGrid& tcv_grid() {
    static const MachineGrid mg = make_machine_grid(load_machine(SHAFRANOV_DATA_DIR "/tcv"),
                                                    Grid(65, 129, 0.55, 1.20, -0.80, 0.80));
    return mg;
}

// Polls until a snapshot with at least min_id has been published, or 30 s pass.
std::shared_ptr<const Snapshot> wait_for(const Worker& w, std::uint64_t min_id) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(30);
    while (std::chrono::steady_clock::now() < deadline) {
        const auto s = w.latest();
        if (s && s->id >= min_id && !w.busy()) {
            return s;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    return w.latest();
}

} // namespace

TEST_CASE("targets_for builds the documented isoflux pairs") {
    Request r;
    const ShapeTargets diverted = targets_for(r);
    CHECK(diverted.isoflux.size() == 11);
    CHECK(diverted.xpoints.size() == 1);

    r.diverted = false;
    const ShapeTargets limited = targets_for(r);
    CHECK(limited.isoflux.size() == 11);
    CHECK(limited.xpoints.empty());

    const std::vector<std::string> fixed = {"OH1", "OH2", "G", "T"};
    CHECK(diverted.fixed == fixed);
    CHECK(limited.fixed == fixed);
}

TEST_CASE("a Worker solves the default request") {
    Worker w(tcv_grid());
    w.request(Request{});
    const auto s = wait_for(w, 1);
    REQUIRE(s);
    CHECK(s->error.empty());
    CHECK(s->eq.converged);
    REQUIRE(s->metrics.has_value());
    CHECK(std::abs(s->metrics->kappa - 1.4) < 0.2);
}

TEST_CASE("only the newest of two quick requests is guaranteed to be solved") {
    Worker w(tcv_grid());
    Request first;
    Request second;
    second.Ip = 1.5e5;
    w.request(first);
    w.request(second);
    const auto s = wait_for(w, 1);
    REQUIRE(s);
    CHECK(s->request.Ip == second.Ip);
    CHECK(s->id <= 2);
}

TEST_CASE("fixed currents of the wrong size give an error snapshot") {
    Worker w(tcv_grid());
    Request r;
    r.fixed_currents = std::vector<double>{1.0, 2.0};
    w.request(r);
    const auto s = wait_for(w, 1);
    REQUIRE(s);
    CHECK_FALSE(s->error.empty());
    CHECK_FALSE(s->metrics.has_value());
}
