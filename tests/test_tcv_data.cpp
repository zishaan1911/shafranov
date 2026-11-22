#include "core/machine.h"

#include <doctest/doctest.h>

#include <algorithm>

TEST_CASE("TCV machine description loads and is self-consistent") {
    const shafranov::Machine m = shafranov::load_machine(SHAFRANOV_DATA_DIR "/tcv");
    CHECK(m.coils.size() == 32);
    CHECK(m.vessel.size() == 256);
    CHECK(m.limiter.size() == 512);

    // The vessel is a parallel set of toroidal paths: Hofmann (1996) quotes 49 uOhm.
    double conductance = 0.0;
    for (const auto& f : m.vessel) {
        conductance += 1.0 / f.resistance;
    }
    CHECK(1.0 / conductance == doctest::Approx(55.7e-6).epsilon(0.01));

    // The first wall lies inside the vessel filaments.
    double vRmin = 1e9, vRmax = 0.0, vZmax = 0.0;
    for (const auto& f : m.vessel) {
        vRmin = std::min(vRmin, f.R);
        vRmax = std::max(vRmax, f.R);
        vZmax = std::max(vZmax, f.Z);
    }
    for (const auto& p : m.limiter) {
        CHECK(p.R > vRmin);
        CHECK(p.R < vRmax);
        CHECK(p.Z < vZmax);
        CHECK(p.Z > -vZmax);
    }

    // G is the fast in-vessel vertical control coil: upper and lower windings in anti-series.
    double g_turns = 0.0;
    for (const auto& c : m.coils) {
        if (c.circuit == "G") {
            g_turns += c.turns;
        }
    }
    CHECK(g_turns == 0.0);
}
