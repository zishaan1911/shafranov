#include "core/constants.h"
#include "core/contour.h"

#include <doctest/doctest.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace shafranov;

namespace {

std::size_t at(const Grid& g, int i, int j) { return static_cast<std::size_t>(g.idx(i, j)); }

// Hand-made topology for an analytic field: nodes with psi > psi_bndry are inside.
Topology hand_topology(const std::vector<double>& psi, double psi_bndry) {
    Topology t;
    t.axis = {0.9, 0.0, 0.0};
    t.psi_bndry = psi_bndry;
    t.limited = false;
    t.inside.resize(psi.size());
    for (std::size_t k = 0; k < psi.size(); ++k) {
        t.inside[k] = psi[k] > psi_bndry ? 1 : 0;
    }
    return t;
}

Grid shape_grid() { return Grid(129, 129, 0.5, 1.3, -0.6, 0.6); }

// Same field as test_topology.cpp: O-point at (0.9, 0), X-point at (0.9, 0.5).
Grid xpoint_grid() { return Grid(65, 129, 0.5, 1.3, -0.8, 0.8); }

std::vector<double> xpoint_psi(const Grid& g) {
    std::vector<double> psi(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        const double x = g.R(i) - 0.9;
        for (int j = 0; j < g.nZ; ++j) {
            const double Z = g.Z(j);
            const double f = -Z * Z / 2.0 + Z * Z * Z / (3.0 * 0.5);
            psi[at(g, i, j)] = -x * x / 2.0 + f;
        }
    }
    return psi;
}

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

std::vector<double> circle_field(const Grid& g) {
    std::vector<double> f(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            const double x = g.R(i) - 1.0;
            f[at(g, i, j)] = x * x + g.Z(j) * g.Z(j);
        }
    }
    return f;
}

} // namespace

TEST_CASE("contour of a circle lies on the circle and has the right length") {
    const Grid g(65, 65, 0.5, 1.5, -0.5, 0.5);
    const std::vector<Segment> segs = contour(g, circle_field(g), 0.09);
    REQUIRE(!segs.empty());
    double length = 0.0;
    for (const Segment& s : segs) {
        CHECK(std::abs(std::hypot(s.R1 - 1.0, s.Z1) - 0.3) < 1e-3);
        CHECK(std::abs(std::hypot(s.R2 - 1.0, s.Z2) - 0.3) < 1e-3);
        length += std::hypot(s.R2 - s.R1, s.Z2 - s.Z1);
    }
    CHECK(length == doctest::Approx(2.0 * pi * 0.3).epsilon(0.01));
}

TEST_CASE("contour: no segments outside the field range, throws on a wrong size") {
    const Grid g(65, 65, 0.5, 1.5, -0.5, 0.5);
    CHECK(contour(g, circle_field(g), 10.0).empty());
    CHECK(contour(g, circle_field(g), -1.0).empty());
    CHECK_THROWS_AS(contour(g, std::vector<double>(10, 0.0), 0.5), std::invalid_argument);
}

TEST_CASE("contour resolves a saddle cell into two segments that do not cross") {
    const Grid g(3, 3, 1.0, 3.0, 0.0, 2.0); // R must be positive; the cell spans R in [1, 2]
    std::vector<double> f(static_cast<std::size_t>(g.size()), 2.0);
    // cell (0, 0) corners alternate 2, 0, 2, 0 with a centre average above the level
    f[at(g, 1, 0)] = 0.0;
    f[at(g, 0, 1)] = 0.0;
    std::vector<Segment> inside_cell;
    for (const Segment& s : contour(g, f, 0.5)) {
        if ((s.R1 + s.R2) / 2.0 < 2.0 && (s.Z1 + s.Z2) / 2.0 < 1.0) {
            inside_cell.push_back(s);
        }
    }
    REQUIRE(inside_cell.size() == 2);
    // the centre is above the level, so the cut-off corners (1, 0) and (0, 1) are separated by
    // the diagonal Z = R - 1
    const bool first_below = inside_cell[0].Z1 < inside_cell[0].R1 - 1.0;
    for (const Segment& s : inside_cell) {
        const bool below = s.Z1 < s.R1 - 1.0;
        CHECK((s.Z2 < s.R2 - 1.0) == below);
    }
    CHECK((inside_cell[1].Z1 < inside_cell[1].R1 - 1.0) != first_below);
}

TEST_CASE("shape metrics of an ellipse") {
    const Grid g = shape_grid();
    std::vector<double> psi(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        const double x = g.R(i) - 0.9;
        for (int j = 0; j < g.nZ; ++j) {
            const double Z = g.Z(j);
            psi[at(g, i, j)] = -(x * x / (0.2 * 0.2) + Z * Z / (0.3 * 0.3));
        }
    }
    const ShapeMetrics m = shape_metrics(g, psi, hand_topology(psi, -1.0));
    CHECK(m.Rgeo == doctest::Approx(0.9).epsilon(1e-2));
    CHECK(m.a == doctest::Approx(0.2).epsilon(1e-2));
    CHECK(m.kappa == doctest::Approx(1.5).epsilon(1e-2));
    CHECK(std::abs(m.delta_upper) < 0.01);
    CHECK(std::abs(m.delta_lower) < 0.01);
}

TEST_CASE("shape metrics: boundary bent toward small R gives positive triangularity") {
    const Grid g = shape_grid();
    std::vector<double> psi(static_cast<std::size_t>(g.size()));
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            const double Z = g.Z(j);
            const double x = (g.R(i) - 0.9) + 0.3 * Z * Z / 0.2;
            psi[at(g, i, j)] = -(x * x / 0.04 + Z * Z / 0.09);
        }
    }
    const ShapeMetrics m = shape_metrics(g, psi, hand_topology(psi, -1.0));
    CHECK(m.delta_upper > 0.0);
    CHECK(m.delta_lower > 0.0);
}

TEST_CASE("shape metrics leave out the private flux region and the legs") {
    const Grid g = xpoint_grid();
    const std::vector<double> psi = xpoint_psi(g);
    const Topology t = find_topology(g, psi, box(0.55, 1.25, -0.75, 0.75));
    const ShapeMetrics m = shape_metrics(g, psi, t);
    CHECK(std::abs(m.Zmax - 0.5) < 2.0 * g.dZ());
    CHECK(std::abs(m.Zmin - (-0.25)) < 0.02);
}
