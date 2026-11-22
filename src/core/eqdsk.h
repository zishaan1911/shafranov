#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace shafranov {

struct GEqdsk {
    std::string description;
    int nw = 0, nh = 0;
    double rdim = 0, zdim = 0, rcentr = 0, rleft = 0, zmid = 0;
    double rmaxis = 0, zmaxis = 0, simag = 0, sibry = 0, bcentr = 0, current = 0;
    std::vector<double> fpol, pres, ffprim, pprime, qpsi; // size nw
    // psirz[i * nh + j]: i over R (0..nw-1), j over Z (0..nh-1). Same layout as Grid fields.
    std::vector<double> psirz;
    std::vector<double> rbbbs, zbbbs, rlim, zlim;
};

// Throws std::invalid_argument if any array size is inconsistent with nw, nh.
void write_geqdsk(std::ostream& out, const GEqdsk& eq);

// Throws std::runtime_error with a clear message on bad input.
GEqdsk read_geqdsk(std::istream& in);

} // namespace shafranov
