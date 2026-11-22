#pragma once

#include "core/grid.h"
#include "core/machine.h"

#include <vector>

namespace shafranov {

struct CriticalPoint {
    double R, Z, psi;
};

struct Topology {
    CriticalPoint axis; // magnetic axis
    // every X-point inside the limiter, nearest in flux to the axis first
    std::vector<CriticalPoint> xpoints;
    double psi_bndry;         // flux on the last closed flux surface
    bool limited;             // true if the limiter, not an X-point, sets psi_bndry
    std::vector<char> inside; // size g.size(): 1 at grid points inside the LCFS
};

// Bilinear interpolation of a grid field at (R, Z); throws std::out_of_range outside the grid.
double interpolate(const Grid& g, const std::vector<double>& f, double R, double Z);

// Point-in-polygon (even-odd rule) for the limiter contour (closed implicitly, last point joins the
// first).
bool inside_polygon(const std::vector<Point>& polygon, double R, double Z);

// Throws std::runtime_error if no O-point lies inside the limiter, std::invalid_argument on a size
// mismatch.
Topology find_topology(const Grid& g, const std::vector<double>& psi,
                       const std::vector<Point>& limiter);

} // namespace shafranov
