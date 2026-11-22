#include "core/profiles.h"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "core/constants.h"
#include "core/grid.h"

namespace {

using shafranov::Grid;
using shafranov::Profile;

constexpr double psi_bndry = -0.04;

struct Setup {
    Grid g{33, 33, 0.6, 1.2, -0.4, 0.4};
    std::vector<double> psi;
    std::vector<char> inside;

    Setup() {
        const auto n = static_cast<std::size_t>(g.size());
        psi.assign(n, 0.0);
        inside.assign(n, 0);
        for (int i = 0; i < g.nR; ++i) {
            for (int j = 0; j < g.nZ; ++j) {
                const double dr = g.R(i) - 0.9;
                const double z = g.Z(j);
                const auto k = static_cast<std::size_t>(g.idx(i, j));
                psi[k] = -(dr * dr + z * z);
                inside[k] = (psi[k] / psi_bndry < 1.0) ? 1 : 0;
            }
        }
    }
};

Profile make_profile(double Ip, double beta0) {
    Profile p;
    p.Ip = Ip;
    p.beta0 = beta0;
    p.R0 = 0.9;
    return p;
}

double shape(const Profile& p, double psi_n) {
    return std::pow(1.0 - std::pow(psi_n, p.alpha_m), p.alpha_n);
}

} // namespace

TEST_CASE("plasma_current integrates to Ip") {
    Setup s;
    for (double Ip : {2e5, -2e5}) {
        const Profile p = make_profile(Ip, 0.5);
        std::vector<double> jphi;
        shafranov::plasma_current(s.g, s.psi, 0.0, psi_bndry, s.inside, p, jphi);
        double sum = 0.0;
        for (double v : jphi) {
            sum += v;
        }
        CHECK(sum * s.g.dR() * s.g.dZ() == doctest::Approx(Ip).epsilon(1e-12));
    }
}

TEST_CASE("jphi is zero outside the mask and positive on the axis") {
    Setup s;
    const Profile p = make_profile(2e5, 0.5);
    std::vector<double> jphi;
    shafranov::plasma_current(s.g, s.psi, 0.0, psi_bndry, s.inside, p, jphi);
    for (std::size_t k = 0; k < jphi.size(); ++k) {
        if (s.inside[k] == 0) {
            CHECK(jphi[k] == 0.0);
        }
    }
    CHECK(jphi[static_cast<std::size_t>(s.g.idx(16, 16))] > 0.0);
}

TEST_CASE("jphi matches the analytic form") {
    Setup s;
    const int i = 18;
    const int j = 17;
    const auto k = static_cast<std::size_t>(s.g.idx(i, j));
    const double R = s.g.R(i);
    const double psi_n = s.psi[k] / psi_bndry;
    for (double beta0 : {0.0, 0.5, 1.0}) {
        const Profile p = make_profile(2e5, beta0);
        std::vector<double> jphi;
        const double lambda =
            shafranov::plasma_current(s.g, s.psi, 0.0, psi_bndry, s.inside, p, jphi);
        const double expected =
            lambda * (beta0 * R / p.R0 + (1.0 - beta0) * p.R0 / R) * shape(p, psi_n);
        CHECK(jphi[k] == doctest::Approx(expected).epsilon(1e-12));
    }
}

TEST_CASE("R p' + FF' / (mu0 R) reproduces jphi") {
    Setup s;
    const Profile p = make_profile(2e5, 0.3);
    std::vector<double> jphi;
    const double lambda = shafranov::plasma_current(s.g, s.psi, 0.0, psi_bndry, s.inside, p, jphi);
    const int pts[5][2] = {{16, 16}, {18, 17}, {14, 15}, {16, 19}, {17, 12}};
    for (const auto& pt : pts) {
        const auto k = static_cast<std::size_t>(s.g.idx(pt[0], pt[1]));
        const double R = s.g.R(pt[0]);
        const double psi_n = s.psi[k] / psi_bndry;
        const double j = R * shafranov::pprime(p, lambda, psi_n) +
                         shafranov::ffprime(p, lambda, psi_n) / (shafranov::mu0 * R);
        CHECK(j == doctest::Approx(jphi[k]).epsilon(1e-12));
    }
}

TEST_CASE("shape is 1 on the axis and 0 at the boundary") {
    const Profile p = make_profile(1.0, 1.0);
    CHECK(shafranov::pprime(p, p.R0, 0.0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(shafranov::pprime(p, p.R0, 1.0) == doctest::Approx(0.0));
}

TEST_CASE("plasma_current rejects bad input") {
    Setup s;
    const Profile p = make_profile(2e5, 0.5);
    std::vector<double> jphi;

    const std::vector<char> none(static_cast<std::size_t>(s.g.size()), 0);
    CHECK_THROWS_AS(shafranov::plasma_current(s.g, s.psi, 0.0, psi_bndry, none, p, jphi),
                    std::runtime_error);
    CHECK_THROWS_AS(shafranov::plasma_current(s.g, s.psi, 0.0, 0.0, s.inside, p, jphi),
                    std::invalid_argument);
    const std::vector<char> short_mask(3, 1);
    CHECK_THROWS_AS(shafranov::plasma_current(s.g, s.psi, 0.0, psi_bndry, short_mask, p, jphi),
                    std::invalid_argument);
}
