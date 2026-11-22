#include "core/constants.h"
#include "core/dst.h"
#include "core/grid.h"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace {

using shafranov::Grid;

std::vector<double> pseudo_random(std::size_t n) {
    std::vector<double> v(n);
    std::uint32_t s = 12345u;
    for (double& e : v) {
        s = s * 1664525u + 1013904223u;
        e = static_cast<double>(s >> 8) / 8388608.0 - 1.0;
    }
    return v;
}

std::vector<double> sample(const Grid& g, double (*f)(double, double)) {
    std::vector<double> psi(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            psi[static_cast<std::size_t>(g.idx(i, j))] = f(g.R(i), g.Z(j));
        }
    }
    return psi;
}

double max_error_r4(int n) {
    const Grid g(n, n, 0.5, 1.5, -1.0, 1.0);
    const std::vector<double> psi = sample(g, [](double R, double) { return R * R * R * R; });
    std::vector<double> out;
    shafranov::delta_star(g, psi, out);
    double err = 0.0;
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            const double R = g.R(i);
            err = std::max(err, std::abs(out[static_cast<std::size_t>(g.idx(i, j))] - 8.0 * R * R));
        }
    }
    return err;
}

} // namespace

TEST_CASE("Grid accessors and validation") {
    const Grid g(5, 3, 1.0, 3.0, -1.0, 1.0);
    CHECK(g.dR() == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(g.dZ() == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(g.R(0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(g.R(4) == doctest::Approx(3.0).epsilon(1e-12));
    CHECK(g.Z(1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(g.idx(2, 1) == 7);
    CHECK(g.size() == 15);

    CHECK_THROWS_AS(Grid(2, 5, 1.0, 2.0, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(Grid(5, 2, 1.0, 2.0, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(Grid(5, 5, 0.0, 2.0, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(Grid(5, 5, 2.0, 2.0, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(Grid(5, 5, 1.0, 2.0, 1.0, 1.0), std::invalid_argument);

    std::vector<double> out;
    CHECK_THROWS_AS(shafranov::delta_star(g, std::vector<double>(3), out), std::invalid_argument);
}

TEST_CASE("delta_star of R^2 vanishes") {
    const Grid g(17, 17, 0.5, 1.5, -1.0, 1.0);
    const std::vector<double> psi = sample(g, [](double R, double) { return R * R; });
    std::vector<double> out;
    shafranov::delta_star(g, psi, out);
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            CHECK(std::abs(out[static_cast<std::size_t>(g.idx(i, j))]) < 1e-9);
        }
    }
}

TEST_CASE("delta_star of Z^2 and R^4") {
    const Grid g(17, 17, 0.5, 1.5, -1.0, 1.0);
    const std::vector<double> psi = sample(g, [](double, double Z) { return Z * Z; });
    std::vector<double> out;
    shafranov::delta_star(g, psi, out);
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            CHECK(out[static_cast<std::size_t>(g.idx(i, j))] == doctest::Approx(2.0).epsilon(1e-9));
        }
    }
    CHECK(out[static_cast<std::size_t>(g.idx(0, 5))] == 0.0);

    CHECK(max_error_r4(17) / max_error_r4(33) > 3.5);
}

TEST_CASE("dst1 matches direct sum") {
    for (std::size_t n : {1u, 3u, 7u, 63u}) {
        const std::vector<double> x = pseudo_random(n);
        std::vector<double> expected(n, 0.0);
        double max_abs = 0.0;
        for (std::size_t k = 1; k <= n; ++k) {
            for (std::size_t j = 1; j <= n; ++j) {
                expected[k - 1] +=
                    x[j - 1] * std::sin(shafranov::pi * static_cast<double>(j * k) /
                                        static_cast<double>(n + 1));
            }
            max_abs = std::max(max_abs, std::abs(expected[k - 1]));
        }
        std::vector<double> got = x;
        shafranov::dst1(got);
        for (std::size_t k = 0; k < n; ++k) {
            CHECK(std::abs(got[k] - expected[k]) < 1e-12 * max_abs);
        }
    }
}

TEST_CASE("dst1 twice scales by (n + 1) / 2") {
    const std::size_t n = 127;
    const std::vector<double> x = pseudo_random(n);
    std::vector<double> y = x;
    shafranov::dst1(y);
    shafranov::dst1(y);
    for (std::size_t k = 0; k < n; ++k) {
        CHECK(y[k] == doctest::Approx(x[k] * static_cast<double>(n + 1) / 2.0).epsilon(1e-9));
    }
}

TEST_CASE("dst1 rejects n + 1 that is not a power of two") {
    std::vector<double> x(6, 1.0);
    CHECK_THROWS_AS(shafranov::dst1(x), std::invalid_argument);
}
