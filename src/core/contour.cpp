#include "core/contour.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace shafranov {

namespace {

struct CellSegment {
    Segment s;
    int i, j; // the cell (i, j)-(i + 1, j + 1) the segment lies in
};

std::vector<CellSegment> march(const Grid& g, const std::vector<double>& f, double level) {
    if (f.size() != static_cast<std::size_t>(g.size())) {
        throw std::invalid_argument("contour: field size does not match the grid");
    }
    std::vector<CellSegment> out;
    for (int i = 0; i + 1 < g.nR; ++i) {
        for (int j = 0; j + 1 < g.nZ; ++j) {
            // corners counter-clockwise from (i, j); edge k joins corner k and corner (k + 1) % 4
            const int ci[4] = {i, i + 1, i + 1, i};
            const int cj[4] = {j, j, j + 1, j + 1};
            double v[4];
            bool above[4];
            double sum = 0.0;
            for (int k = 0; k < 4; ++k) {
                v[k] = f[static_cast<std::size_t>(g.idx(ci[k], cj[k]))];
                above[k] = v[k] >= level;
                sum += v[k];
            }

            const auto emit = [&](int ea, int eb) {
                double p[2][2];
                const int edges[2] = {ea, eb};
                for (int m = 0; m < 2; ++m) {
                    const int a = edges[m];
                    const int b = (a + 1) % 4;
                    const double t = (level - v[a]) / (v[b] - v[a]);
                    const double Ra = g.R(ci[a]);
                    const double Za = g.Z(cj[a]);
                    p[m][0] = Ra + t * (g.R(ci[b]) - Ra);
                    p[m][1] = Za + t * (g.Z(cj[b]) - Za);
                }
                out.push_back({{p[0][0], p[0][1], p[1][0], p[1][1]}, i, j});
            };

            int crossed[4];
            int n = 0;
            for (int e = 0; e < 4; ++e) {
                if (above[e] != above[(e + 1) % 4]) {
                    crossed[n++] = e;
                }
            }
            if (n == 2) {
                emit(crossed[0], crossed[1]);
            } else if (n == 4) {
                // saddle: corners that disagree with the cell centre are cut off one by one
                const bool centre = sum / 4.0 >= level;
                for (int k = 0; k < 4; ++k) {
                    if (above[k] != centre) {
                        emit((k + 3) % 4, k);
                    }
                }
            }
        }
    }
    return out;
}

} // namespace

std::vector<Segment> contour(const Grid& g, const std::vector<double>& f, double level) {
    std::vector<Segment> segs;
    for (const CellSegment& cs : march(g, f, level)) {
        segs.push_back(cs.s);
    }
    return segs;
}

ShapeMetrics shape_metrics(const Grid& g, const std::vector<double>& psi, const Topology& t) {
    bool any = false;
    ShapeMetrics m{};
    for (const CellSegment& cs : march(g, psi, t.psi_bndry)) {
        const int ci[4] = {cs.i, cs.i + 1, cs.i + 1, cs.i};
        const int cj[4] = {cs.j, cs.j, cs.j + 1, cs.j + 1};
        bool touches = false;
        for (int k = 0; k < 4; ++k) {
            touches = touches || t.inside[static_cast<std::size_t>(g.idx(ci[k], cj[k]))] != 0;
        }
        if (!touches) {
            continue;
        }
        const double Rs[2] = {cs.s.R1, cs.s.R2};
        const double Zs[2] = {cs.s.Z1, cs.s.Z2};
        for (int k = 0; k < 2; ++k) {
            if (!any) {
                m.Rmin = m.Rmax = Rs[k];
                m.Zmin = m.Zmax = Zs[k];
                m.R_at_Zmax = m.R_at_Zmin = Rs[k];
                any = true;
                continue;
            }
            if (Rs[k] < m.Rmin) {
                m.Rmin = Rs[k];
            }
            if (Rs[k] > m.Rmax) {
                m.Rmax = Rs[k];
            }
            if (Zs[k] > m.Zmax) {
                m.Zmax = Zs[k];
                m.R_at_Zmax = Rs[k];
            }
            if (Zs[k] < m.Zmin) {
                m.Zmin = Zs[k];
                m.R_at_Zmin = Rs[k];
            }
        }
    }
    if (!any) {
        throw std::runtime_error("shape_metrics: no boundary segment touches the plasma");
    }
    m.Rgeo = (m.Rmax + m.Rmin) / 2.0;
    m.a = (m.Rmax - m.Rmin) / 2.0;
    m.kappa = (m.Zmax - m.Zmin) / (2.0 * m.a);
    m.delta_upper = (m.Rgeo - m.R_at_Zmax) / m.a;
    m.delta_lower = (m.Rgeo - m.R_at_Zmin) / m.a;
    return m;
}

} // namespace shafranov
