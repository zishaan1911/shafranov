#include "core/coils.h"

#include <doctest/doctest.h>

#include <cmath>

TEST_CASE("TCV circuits keep their turns when split into filaments") {
    const shafranov::Machine m = shafranov::load_machine(SHAFRANOV_DATA_DIR "/tcv");
    const auto names = shafranov::circuit_names(m);
    const auto windings = shafranov::circuit_windings(m);
    REQUIRE(names.size() == 20);
    REQUIRE(windings.size() == names.size());
    CHECK(names.front() == "OH1");

    double oh1 = 0.0, g = 0.0;
    for (const auto& w : windings[0]) {
        oh1 += w.turns;
        CHECK(std::abs(w.R - 0.4225) <= 0.0315);
    }
    for (std::size_t k = 0; k < names.size(); ++k) {
        if (names[k] == "G") {
            for (const auto& w : windings[k]) {
                g += w.turns;
            }
        }
    }
    CHECK(oh1 == doctest::Approx(143.0));
    CHECK(windings[0].size() == 4 * 80); // 0.063 x 1.584 m in 2 cm pieces
    CHECK(g == doctest::Approx(0.0));
}
