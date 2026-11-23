#include "core/topology.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace shafranov {

namespace {

struct Candidate {
    double R, Z, psi;
    bool xpoint;
    bool maximum; // meaningful for O-points: psi has a local maximum
};

struct Derivs {
    double pR, pZ, pRR, pZZ, pRZ;
};

std::size_t at(const Grid& g, int i, int j) { return static_cast<std::size_t>(g.idx(i, j)); }

Derivs derivs(const Grid& g, const std::vector<double>& psi, int i, int j) {
    const double dR = g.dR();
    const double dZ = g.dZ();
    const double c = psi[at(g, i, j)];
    const double rp = psi[at(g, i + 1, j)];
    const double rm = psi[at(g, i - 1, j)];
    const double zp = psi[at(g, i, j + 1)];
    const double zm = psi[at(g, i, j - 1)];
    Derivs d;
    d.pR = (rp - rm) / (2.0 * dR);
    d.pZ = (zp - zm) / (2.0 * dZ);
    d.pRR = (rp - 2.0 * c + rm) / (dR * dR);
    d.pZZ = (zp - 2.0 * c + zm) / (dZ * dZ);
    d.pRZ = (psi[at(g, i + 1, j + 1)] - psi[at(g, i + 1, j - 1)] - psi[at(g, i - 1, j + 1)] +
             psi[at(g, i - 1, j - 1)]) /
            (4.0 * dR * dZ);
    return d;
}

int nearest_index(double v, double vmin, double d, int n) {
    const int k = static_cast<int>(std::lround((v - vmin) / d));
    return std::clamp(k, 0, n - 1);
}

bool in_grid(const Grid& g, double R, double Z) {
    return R >= g.Rmin && R <= g.Rmax && Z >= g.Zmin && Z <= g.Zmax;
}

// Grid points connected to the axis with psi_n < 1, cut off beyond the first X-point if given.
std::vector<char> region(const Grid& g, const std::vector<double>& psi, const CriticalPoint& axis,
                         const CriticalPoint* x1, double psi_b) {
    std::vector<char> in(static_cast<std::size_t>(g.size()), 0);

    auto ok = [&](int i, int j) {
        const double pn = (psi[at(g, i, j)] - axis.psi) / (psi_b - axis.psi);
        if (!(pn < 1.0)) {
            return false;
        }
        if (x1 != nullptr) {
            const double dot = (g.R(i) - x1->R) * (axis.R - x1->R) +
                               (g.Z(j) - x1->Z) * (axis.Z - x1->Z);
            if (dot < 0.0) {
                return false;
            }
        }
        return true;
    };

    const int si = nearest_index(axis.R, g.Rmin, g.dR(), g.nR);
    const int sj = nearest_index(axis.Z, g.Zmin, g.dZ(), g.nZ);
    if (!ok(si, sj)) {
        return in;
    }
    in[at(g, si, sj)] = 1;
    std::vector<std::pair<int, int>> stack;
    stack.emplace_back(si, sj);
    const int di[4] = {1, -1, 0, 0};
    const int dj[4] = {0, 0, 1, -1};
    while (!stack.empty()) {
        const auto [i, j] = stack.back();
        stack.pop_back();
        for (int n = 0; n < 4; ++n) {
            const int ni = i + di[n];
            const int nj = j + dj[n];
            if (ni < 0 || ni >= g.nR || nj < 0 || nj >= g.nZ) {
                continue;
            }
            if (in[at(g, ni, nj)] != 0 || !ok(ni, nj)) {
                continue;
            }
            in[at(g, ni, nj)] = 1;
            stack.emplace_back(ni, nj);
        }
    }
    return in;
}

} // namespace

double interpolate(const Grid& g, const std::vector<double>& f, double R, double Z) {
    if (!in_grid(g, R, Z)) {
        throw std::out_of_range("interpolate: point outside the grid");
    }
    const double x = (R - g.Rmin) / g.dR();
    const double y = (Z - g.Zmin) / g.dZ();
    const int i = std::min(static_cast<int>(x), g.nR - 2);
    const int j = std::min(static_cast<int>(y), g.nZ - 2);
    const double a = x - i;
    const double b = y - j;
    return (1.0 - a) * (1.0 - b) * f[at(g, i, j)] + a * (1.0 - b) * f[at(g, i + 1, j)] +
           (1.0 - a) * b * f[at(g, i, j + 1)] + a * b * f[at(g, i + 1, j + 1)];
}

bool inside_polygon(const std::vector<Point>& polygon, double R, double Z) {
    bool in = false;
    const std::size_t n = polygon.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const Point& a = polygon[i];
        const Point& b = polygon[j];
        if ((a.Z > Z) != (b.Z > Z) && R < (b.R - a.R) * (Z - a.Z) / (b.Z - a.Z) + a.R) {
            in = !in;
        }
    }
    return in;
}

Topology find_topology(const Grid& g, const std::vector<double>& psi,
                       const std::vector<Point>& limiter, const Point* axis_hint) {
    if (psi.size() != static_cast<std::size_t>(g.size())) {
        throw std::invalid_argument("find_topology: psi size does not match grid");
    }
    const double dR = g.dR();
    const double dZ = g.dZ();

    // |grad psi|^2 at interior nodes (zero elsewhere, never read there).
    std::vector<double> grad2(static_cast<std::size_t>(g.size()), 0.0);
    for (int i = 1; i < g.nR - 1; ++i) {
        for (int j = 1; j < g.nZ - 1; ++j) {
            const Derivs d = derivs(g, psi, i, j);
            grad2[at(g, i, j)] = d.pR * d.pR + d.pZ * d.pZ;
        }
    }

    std::vector<Candidate> found;
    for (int i = 2; i < g.nR - 2; ++i) {
        for (int j = 2; j < g.nZ - 2; ++j) {
            const double c = grad2[at(g, i, j)];
            bool smallest = true;
            for (int a = -1; a <= 1 && smallest; ++a) {
                for (int b = -1; b <= 1; ++b) {
                    if ((a != 0 || b != 0) && grad2[at(g, i + a, j + b)] < c) {
                        smallest = false;
                        break;
                    }
                }
            }
            if (!smallest) {
                continue;
            }

            const Derivs d = derivs(g, psi, i, j);
            const double det = d.pRR * d.pZZ - d.pRZ * d.pRZ;
            if (det == 0.0) {
                continue;
            }
            const double deltaR = (-d.pR * d.pZZ + d.pZ * d.pRZ) / det;
            const double deltaZ = (-d.pZ * d.pRR + d.pR * d.pRZ) / det;
            if (std::abs(deltaR) > dR || std::abs(deltaZ) > dZ) {
                continue;
            }
            Candidate p;
            p.R = g.R(i) + deltaR;
            p.Z = g.Z(j) + deltaZ;
            p.psi = psi[at(g, i, j)] + d.pR * deltaR + d.pZ * deltaZ +
                    0.5 * (d.pRR * deltaR * deltaR + 2.0 * d.pRZ * deltaR * deltaZ +
                           d.pZZ * deltaZ * deltaZ);
            p.xpoint = det < 0.0;
            p.maximum = d.pRR + d.pZZ < 0.0;

            bool duplicate = false;
            for (const Candidate& q : found) {
                if (std::abs(q.R - p.R) < dR && std::abs(q.Z - p.Z) < dZ) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                found.push_back(p);
            }
        }
    }

    const double Rc = axis_hint != nullptr ? axis_hint->R : 0.5 * (g.Rmin + g.Rmax);
    const double Zc = axis_hint != nullptr ? axis_hint->Z : 0.5 * (g.Zmin + g.Zmax);
    const Candidate* best = nullptr;
    double best_dist = 0.0;
    for (const Candidate& p : found) {
        if (p.xpoint || !inside_polygon(limiter, p.R, p.Z)) {
            continue;
        }
        const double dist = (p.R - Rc) * (p.R - Rc) + (p.Z - Zc) * (p.Z - Zc);
        if (best == nullptr || dist < best_dist) {
            best = &p;
            best_dist = dist;
        }
    }
    if (best == nullptr) {
        throw std::runtime_error("find_topology: no magnetic axis inside the limiter");
    }

    Topology t;
    t.axis = {best->R, best->Z, best->psi};
    const bool axis_max = best->maximum;

    for (const Candidate& p : found) {
        if (p.xpoint && inside_polygon(limiter, p.R, p.Z)) {
            t.xpoints.push_back({p.R, p.Z, p.psi});
        }
    }
    std::sort(t.xpoints.begin(), t.xpoints.end(),
              [&](const CriticalPoint& a, const CriticalPoint& b) {
                  return std::abs(a.psi - t.axis.psi) < std::abs(b.psi - t.axis.psi);
              });

    // Keeps the limiter flux nearest the axis: the largest if the axis is a maximum, else the
    // smallest.
    bool have_flux = false;
    double limiter_flux = 0.0;
    auto consider = [&](double v) {
        if (!have_flux || (axis_max ? v > limiter_flux : v < limiter_flux)) {
            limiter_flux = v;
            have_flux = true;
        }
    };

    const CriticalPoint* x1 = t.xpoints.empty() ? nullptr : &t.xpoints.front();
    if (x1 != nullptr) {
        const std::vector<char> core = region(g, psi, t.axis, x1, x1->psi);
        for (const Point& p : limiter) {
            if (!in_grid(g, p.R, p.Z)) {
                continue;
            }
            const int i = nearest_index(p.R, g.Rmin, dR, g.nR);
            const int j = nearest_index(p.Z, g.Zmin, dZ, g.nZ);
            if (core[at(g, i, j)] != 0) {
                consider(interpolate(g, psi, p.R, p.Z));
            }
        }
        t.limited = have_flux;
        t.psi_bndry = have_flux ? limiter_flux : x1->psi;
    } else {
        for (const Point& p : limiter) {
            if (in_grid(g, p.R, p.Z)) {
                consider(interpolate(g, psi, p.R, p.Z));
            }
        }
        if (!have_flux) {
            throw std::runtime_error("find_topology: limiter lies outside the grid");
        }
        t.limited = true;
        t.psi_bndry = limiter_flux;
    }

    t.inside = region(g, psi, t.axis, x1, t.psi_bndry);
    for (int i = 0; i < g.nR; ++i) {
        for (int j = 0; j < g.nZ; ++j) {
            char& flag = t.inside[at(g, i, j)];
            if (flag != 0 && !inside_polygon(limiter, g.R(i), g.Z(j))) {
                flag = 0;
            }
        }
    }
    return t;
}

} // namespace shafranov
