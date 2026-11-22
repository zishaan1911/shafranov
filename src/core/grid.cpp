#include "core/grid.h"

#include <cstddef>
#include <stdexcept>

namespace shafranov {

Grid::Grid(int nR_, int nZ_, double Rmin_, double Rmax_, double Zmin_, double Zmax_)
    : nR(nR_), nZ(nZ_), Rmin(Rmin_), Rmax(Rmax_), Zmin(Zmin_), Zmax(Zmax_) {
    if (nR < 3 || nZ < 3) {
        throw std::invalid_argument("Grid: nR and nZ must be at least 3");
    }
    if (Rmin <= 0.0) {
        throw std::invalid_argument("Grid: Rmin must be positive");
    }
    if (Rmax <= Rmin || Zmax <= Zmin) {
        throw std::invalid_argument("Grid: require Rmax > Rmin and Zmax > Zmin");
    }
}

double Grid::dR() const { return (Rmax - Rmin) / (nR - 1); }

double Grid::dZ() const { return (Zmax - Zmin) / (nZ - 1); }

double Grid::R(int i) const { return Rmin + i * dR(); }

double Grid::Z(int j) const { return Zmin + j * dZ(); }

void delta_star(const Grid& g, const std::vector<double>& psi, std::vector<double>& out) {
    if (psi.size() != static_cast<std::size_t>(g.size())) {
        throw std::invalid_argument("delta_star: psi size does not match grid");
    }
    out.assign(static_cast<std::size_t>(g.size()), 0.0);

    const double dR = g.dR();
    const double dZ = g.dZ();
    for (int i = 1; i < g.nR - 1; ++i) {
        const double R = g.R(i);
        for (int j = 1; j < g.nZ - 1; ++j) {
            const double c = psi[static_cast<std::size_t>(g.idx(i, j))];
            const double rp = psi[static_cast<std::size_t>(g.idx(i + 1, j))];
            const double rm = psi[static_cast<std::size_t>(g.idx(i - 1, j))];
            const double zp = psi[static_cast<std::size_t>(g.idx(i, j + 1))];
            const double zm = psi[static_cast<std::size_t>(g.idx(i, j - 1))];
            out[static_cast<std::size_t>(g.idx(i, j))] =
                (rp - 2.0 * c + rm) / (dR * dR) - (rp - rm) / (2.0 * R * dR) +
                (zp - 2.0 * c + zm) / (dZ * dZ);
        }
    }
}

} // namespace shafranov
