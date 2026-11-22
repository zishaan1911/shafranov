#include "core/topology.h"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace shafranov;

namespace {

Grid test_grid() { return Grid(65, 129, 0.5, 1.3, -0.8, 0.8); }

// O-point (maximum) at (0.9, 0), X-point at (0.9, 0.5); sign = -1 flips the field.
std::vector<double> test_psi(const Grid& g, double sign) {
    std::vector<double> psi(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        const double x = g.R(i) - 0.9;
        for (int j = 0; j < g.nZ; ++j) {
            const double Z = g.Z(j);
            const double f = -Z * Z / 2.0 + Z * Z * Z / (3.0 * 0.5);
            psi[static_cast<std::size_t>(g.idx(i, j))] = sign * (-x * x / 2.0 + f);
        }
    }
    return psi;
}

// Rectangle sampled with 40 points per edge, in order around the box.
std::vector<Point> box(double Rlo, double Rhi, double Zlo, double Zhi) {
    const int n = 40;
    std::vector<Point> pts;
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rlo + (Rhi - Rlo) * k / n, Zlo});
    }
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rhi, Zlo + (Zhi - Zlo) * k / n});
    }
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rhi - (Rhi - Rlo) * k / n, Zhi});
    }
    for (int k = 0; k < n; ++k) {
        pts.push_back({Rlo, Zhi - (Zhi - Zlo) * k / n});
    }
    return pts;
}

std::size_t node(const Grid& g, double R, double Z) {
    const int i = static_cast<int>(std::lround((R - g.Rmin) / g.dR()));
    const int j = static_cast<int>(std::lround((Z - g.Zmin) / g.dZ()));
    return static_cast<std::size_t>(g.idx(i, j));
}

void check_diverted(double sign) {
    const Grid g = test_grid();
    const Topology t = find_topology(g, test_psi(g, sign), box(0.55, 1.25, -0.75, 0.75));
    CHECK(std::abs(t.axis.R - 0.9) < 0.1 * g.dR());
    CHECK(std::abs(t.axis.Z - 0.0) < 0.1 * g.dR());
    CHECK(std::abs(t.axis.psi) < 1e-6);
    REQUIRE(t.xpoints.size() == 1);
    CHECK(std::abs(t.xpoints[0].R - 0.9) < 0.5 * g.dR());
    CHECK(std::abs(t.xpoints[0].Z - 0.5) < 0.5 * g.dZ());
    CHECK(std::abs(t.psi_bndry - sign * (-1.0 / 24.0)) < 1e-4);
    CHECK(t.limited == false);
    REQUIRE(t.inside.size() == static_cast<std::size_t>(g.size()));
    CHECK(t.inside[node(g, t.axis.R, t.axis.Z)] == 1);
    CHECK(t.inside[node(g, 0.9, 0.6)] == 0);
}

void check_limited(double sign) {
    const Grid g = test_grid();
    const Topology t = find_topology(g, test_psi(g, sign), box(0.7, 1.1, -0.4, 0.4));
    CHECK(t.limited == true);
    CHECK(std::abs(t.psi_bndry - sign * -0.02) < 1e-4);
    CHECK(std::abs(t.axis.R - 0.9) < 0.1 * g.dR());
    CHECK(t.inside[node(g, 0.9, 0.0)] == 1);
}

} // namespace

TEST_CASE("diverted plasma: X-point sets the boundary") { check_diverted(1.0); }

TEST_CASE("limited plasma: limiter sets the boundary") { check_limited(1.0); }

TEST_CASE("flipping the sign of psi gives the same topology") {
    check_diverted(-1.0);
    check_limited(-1.0);
}

TEST_CASE("interpolate is exact for bilinear fields") {
    const Grid g = test_grid();
    std::vector<double> f(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            f[static_cast<std::size_t>(g.idx(i, j))] = 2.0 + 3.0 * g.R(i) - 4.0 * g.Z(j) +
                                                       0.5 * g.R(i) * g.Z(j);
        }
    }
    const double R = 0.8123;
    const double Z = 0.1234;
    CHECK(interpolate(g, f, R, Z) == doctest::Approx(2.0 + 3.0 * R - 4.0 * Z + 0.5 * R * Z));
    CHECK_THROWS_AS(interpolate(g, f, 0.4, 0.0), std::out_of_range);
    CHECK_THROWS_AS(interpolate(g, f, 0.9, 0.9), std::out_of_range);
}

TEST_CASE("inside_polygon handles convex and concave polygons") {
    const std::vector<Point> square = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    CHECK(inside_polygon(square, 0.5, 0.5));
    CHECK_FALSE(inside_polygon(square, 1.5, 0.5));
    CHECK_FALSE(inside_polygon(square, 0.5, -0.5));

    const std::vector<Point> ell = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0},
                                    {1.0, 1.0}, {1.0, 2.0}, {0.0, 2.0}};
    CHECK(inside_polygon(ell, 0.5, 1.5));
    CHECK(inside_polygon(ell, 1.5, 0.5));
    CHECK_FALSE(inside_polygon(ell, 1.5, 1.5));
}

TEST_CASE("find_topology throws when the limiter excludes the O-point") {
    const Grid g = test_grid();
    CHECK_THROWS_AS(find_topology(g, test_psi(g, 1.0), box(1.0, 1.2, -0.1, 0.1)),
                    std::runtime_error);
}
