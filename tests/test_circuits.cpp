#include "core/circuits.h"

#include "core/constants.h"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

namespace {

using namespace shafranov;

struct Fixture {
    Machine machine;
    MachineGrid mg;
    CircuitModel model;
};

const Fixture& tcv() {
    static const Fixture f = [] {
        Machine machine = load_machine(SHAFRANOV_DATA_DIR "/tcv");
        MachineGrid mg = make_machine_grid(machine, Grid(33, 65, 0.55, 1.20, -0.80, 0.80));
        CircuitModel model = make_circuit_model(machine, mg);
        return Fixture{std::move(machine), std::move(mg), std::move(model)};
    }();
    return f;
}

int index_of(const std::string& name) {
    const auto& names = tcv().mg.circuits;
    const auto it = std::find(names.begin(), names.end(), name);
    REQUIRE(it != names.end());
    return static_cast<int>(it - names.begin());
}

} // namespace

TEST_CASE("self_inductance of a thin ring matches the closed form") {
    CHECK(self_inductance(1.0, 0.01, 0.01) ==
          doctest::Approx(mu0 * (std::log(8.0 / (0.2235 * 0.02)) - 2.0)).epsilon(1e-12));
}

TEST_CASE("the TCV inductance matrix is symmetric and positive definite") {
    const CircuitModel& c = tcv().model;
    const Eigen::Index n = c.M.rows();
    REQUIRE(n == c.n_circuits + c.n_vessel);
    REQUIRE(c.M.cols() == n);
    REQUIRE(c.R.size() == n);
    for (Eigen::Index i = 0; i < n; ++i) {
        CHECK(c.M(i, i) > 0.0);
        for (Eigen::Index j = 0; j < i; ++j) {
            CHECK(c.M(i, j) == c.M(j, i));
        }
    }
    const Eigen::LLT<Eigen::MatrixXd> llt(c.M);
    CHECK(llt.info() == Eigen::Success);
}

TEST_CASE("TCV circuit inductances and resistances match Marchioni (2024)") {
    const CircuitModel& c = tcv().model;
    const auto L = [&](const char* name) {
        const int k = index_of(name);
        return c.M(k, k);
    };
    CHECK(L("OH1") == doctest::Approx(6.9e-3).epsilon(0.05));
    CHECK(L("OH2") == doctest::Approx(4.8e-3).epsilon(0.05));
    CHECK(L("E1") == doctest::Approx(1.7e-3).epsilon(0.05));
    CHECK(L("F1") == doctest::Approx(7.4e-3).epsilon(0.05));
    CHECK(L("G") == doctest::Approx(0.1e-3).epsilon(0.10));
    CHECK(c.R(index_of("OH1")) == doctest::Approx(15e-3).epsilon(0.05));
    CHECK(c.R(index_of("F1")) == doctest::Approx(39.1e-3).epsilon(0.05));
}

TEST_CASE("TCV vessel modes match Hofmann (1996)") {
    const Fixture& f = tcv();
    const auto modes = vessel_modes(f.model, f.machine);
    REQUIRE(static_cast<int>(modes.size()) == f.model.n_vessel);

    double longest_even = 0.0, longest_odd = 0.0;
    for (const VesselMode& md : modes) {
        CHECK(std::isfinite(md.tau));
        CHECK(md.tau > 0.0);
        if (md.parity > 0 && longest_even == 0.0) {
            longest_even = md.tau;
        }
        if (md.parity < 0 && longest_odd == 0.0) {
            longest_odd = md.tau;
        }
    }
    CHECK(modes.front().parity == 1);
    MESSAGE("longest even mode: " << longest_even * 1e3 << " ms");
    MESSAGE("longest odd mode: " << longest_odd * 1e3 << " ms");
    CHECK(longest_even == doctest::Approx(14.0e-3).epsilon(0.10));
    CHECK(longest_odd == doctest::Approx(8.6e-3).epsilon(0.10));
}
