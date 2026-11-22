#pragma once

#include <vector>

namespace shafranov {

// Uniform rectangular grid in (R, Z) including boundary points. Fields are stored as
// std::vector<double> of size nR * nZ, index i * nZ + j (i over R, j over Z), so a column of
// constant R is contiguous (the sine transform runs along Z).
struct Grid {
    int nR, nZ;
    double Rmin, Rmax, Zmin, Zmax;

    Grid(int nR, int nZ, double Rmin, double Rmax, double Zmin, double Zmax); // validates

    double dR() const;     // (Rmax - Rmin) / (nR - 1)
    double dZ() const;     // (Zmax - Zmin) / (nZ - 1)
    double R(int i) const; // Rmin + i * dR()
    double Z(int j) const; // Zmin + j * dZ()
    int idx(int i, int j) const { return i * nZ + j; }
    int size() const { return nR * nZ; }
};

// out = Delta* psi at interior points, 0 on the boundary. out is resized to g.size().
// Delta* psi = d2psi/dR2 - (1/R) dpsi/dR + d2psi/dZ2, second-order central differences.
void delta_star(const Grid& g, const std::vector<double>& psi, std::vector<double>& out);

} // namespace shafranov
