#include "core/eqdsk.h"

#include <doctest/doctest.h>

#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using shafranov::GEqdsk;

GEqdsk make_equilibrium() {
    GEqdsk eq;
    eq.description = "round trip test";
    eq.nw = 7;
    eq.nh = 5;
    eq.rdim = 1.7;
    eq.zdim = 3.2;
    eq.rcentr = 1.69;
    eq.rleft = 0.84;
    eq.zmid = -0.05;
    eq.rmaxis = 1.65;
    eq.zmaxis = 0.02;
    eq.simag = -0.4;
    eq.sibry = 0.13;
    eq.bcentr = 2.1;
    eq.current = 1.5e6;
    const std::size_t nw = 7;
    for (std::size_t i = 0; i < nw; ++i) {
        const double x = static_cast<double>(i);
        eq.fpol.push_back(3.5 + 0.1 * x);
        eq.pres.push_back(1.0e4 - 1.3e3 * x);
        eq.ffprim.push_back(-0.7 + 0.05 * x);
        eq.pprime.push_back(-2.0e4 * (1.0 + x));
        eq.qpsi.push_back(1.1 + 0.4 * x);
    }
    for (int i = 0; i < eq.nw; ++i)
        for (int j = 0; j < eq.nh; ++j)
            eq.psirz.push_back(i + 10.0 * j - 20.25);
    eq.rbbbs = {1.2, 1.9, 1.5};
    eq.zbbbs = {0.1, -0.3, 0.8};
    eq.rlim = {0.9, 2.4, 2.4, 0.9};
    eq.zlim = {-1.6, -1.6, 1.6, 1.6};
    return eq;
}

void check_vectors(const std::vector<double>& a, const std::vector<double>& b) {
    REQUIRE(a.size() == b.size());
    for (std::size_t k = 0; k < a.size(); ++k)
        CHECK(a[k] == doctest::Approx(b[k]).epsilon(1e-9));
}

} // namespace

TEST_CASE("G-EQDSK round trip preserves every field") {
    const GEqdsk eq = make_equilibrium();
    std::stringstream ss;
    shafranov::write_geqdsk(ss, eq);
    const GEqdsk back = shafranov::read_geqdsk(ss);

    CHECK(back.description == eq.description);
    CHECK(back.nw == eq.nw);
    CHECK(back.nh == eq.nh);
    CHECK(back.rdim == doctest::Approx(eq.rdim).epsilon(1e-9));
    CHECK(back.zdim == doctest::Approx(eq.zdim).epsilon(1e-9));
    CHECK(back.rcentr == doctest::Approx(eq.rcentr).epsilon(1e-9));
    CHECK(back.rleft == doctest::Approx(eq.rleft).epsilon(1e-9));
    CHECK(back.zmid == doctest::Approx(eq.zmid).epsilon(1e-9));
    CHECK(back.rmaxis == doctest::Approx(eq.rmaxis).epsilon(1e-9));
    CHECK(back.zmaxis == doctest::Approx(eq.zmaxis).epsilon(1e-9));
    CHECK(back.simag == doctest::Approx(eq.simag).epsilon(1e-9));
    CHECK(back.sibry == doctest::Approx(eq.sibry).epsilon(1e-9));
    CHECK(back.bcentr == doctest::Approx(eq.bcentr).epsilon(1e-9));
    CHECK(back.current == doctest::Approx(eq.current).epsilon(1e-9));
    check_vectors(back.fpol, eq.fpol);
    check_vectors(back.pres, eq.pres);
    check_vectors(back.ffprim, eq.ffprim);
    check_vectors(back.pprime, eq.pprime);
    check_vectors(back.qpsi, eq.qpsi);
    check_vectors(back.psirz, eq.psirz);
    check_vectors(back.rbbbs, eq.rbbbs);
    check_vectors(back.zbbbs, eq.zbbbs);
    check_vectors(back.rlim, eq.rlim);
    check_vectors(back.zlim, eq.zlim);
}

TEST_CASE("G-EQDSK header and scalar lines have fixed widths") {
    std::stringstream ss;
    shafranov::write_geqdsk(ss, make_equilibrium());
    std::string line1, line2;
    std::getline(ss, line1);
    std::getline(ss, line2);
    CHECK(line1.size() == 60);
    CHECK(line2.size() == 80);
}

TEST_CASE("G-EQDSK reader parses run-together negative fields") {
    const std::string one = " 1.000000000E+00";
    const std::string five = one + one + one + one + one + "\n";
    const std::string text = std::string("hand-written").append(36, ' ') + "   0   2   2\n" +
                             five + five + five + five +
                             "-2.000000000E+00-3.000000000E+00\n" // fpol
                             + one + one + "\r\n"                 // pres, CRLF line ending
                             + one + one + "\n"                   // ffprim
                             + one + one + "\n"                   // pprime
                             + one + one + one + one + "\n"       // psirz
                             + one + one + "\n"                   // qpsi
                             + "    1    1\n" + one + one + "\n" + one + one + "\n";
    std::istringstream in(text);
    const GEqdsk eq = shafranov::read_geqdsk(in);
    CHECK(eq.description == "hand-written");
    CHECK(eq.nw == 2);
    CHECK(eq.nh == 2);
    REQUIRE(eq.fpol.size() == 2);
    CHECK(eq.fpol[0] == doctest::Approx(-2.0).epsilon(1e-12));
    CHECK(eq.fpol[1] == doctest::Approx(-3.0).epsilon(1e-12));
    CHECK(eq.pres.size() == 2);
    CHECK(eq.psirz.size() == 4);
    CHECK(eq.rbbbs.size() == 1);
    CHECK(eq.zlim.size() == 1);
}

TEST_CASE("G-EQDSK reader throws on a truncated file") {
    std::stringstream ss;
    shafranov::write_geqdsk(ss, make_equilibrium());
    const std::string full = ss.str();
    std::istringstream in(full.substr(0, full.size() / 2));
    CHECK_THROWS_AS(shafranov::read_geqdsk(in), std::runtime_error);
}

TEST_CASE("G-EQDSK writer rejects inconsistent array sizes") {
    GEqdsk eq = make_equilibrium();
    eq.fpol.pop_back();
    std::stringstream ss;
    CHECK_THROWS_AS(shafranov::write_geqdsk(ss, eq), std::invalid_argument);
}
