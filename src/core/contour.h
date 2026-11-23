#pragma once

#include "core/grid.h"
#include "core/topology.h"

#include <vector>

namespace shafranov {

struct Segment {
    double R1, Z1, R2, Z2;
};

// Marching squares: the line segments where the grid field f crosses level, with linear
// interpolation along cell edges. Saddle cells are resolved by the cell-centre average.
std::vector<Segment> contour(const Grid& g, const std::vector<double>& f, double level);

// Geometry of the last closed flux surface, from the psi_bndry contour segments that lie in cells
// touching t.inside (so the separatrix legs and private flux region are left out).
struct ShapeMetrics {
    double Rmin, Rmax, Zmin, Zmax;
    double R_at_Zmax, R_at_Zmin;
    double Rgeo;        // (Rmax + Rmin) / 2
    double a;           // minor radius (Rmax - Rmin) / 2
    double kappa;       // elongation (Zmax - Zmin) / (2 a)
    double delta_upper; // triangularity (Rgeo - R_at_Zmax) / a
    double delta_lower; // (Rgeo - R_at_Zmin) / a
};

// Throws std::runtime_error if no boundary segment touches the plasma.
ShapeMetrics shape_metrics(const Grid& g, const std::vector<double>& psi, const Topology& t);

} // namespace shafranov
