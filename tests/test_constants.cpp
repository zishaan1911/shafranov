#include "core/constants.h"

#include <doctest/doctest.h>

TEST_CASE("mu0 is close to 4 pi 1e-7") {
    CHECK(shafranov::mu0 == doctest::Approx(4e-7 * shafranov::pi).epsilon(1e-9));
}
