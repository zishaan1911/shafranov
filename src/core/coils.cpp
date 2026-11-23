#include "core/coils.h"

#include <algorithm>
#include <cmath>

namespace shafranov {

std::vector<std::string> circuit_names(const Machine& m) {
    std::vector<std::string> names;
    for (const Coil& c : m.coils) {
        if (std::find(names.begin(), names.end(), c.circuit) == names.end()) {
            names.push_back(c.circuit);
        }
    }
    return names;
}

std::vector<std::vector<Winding>> circuit_windings(const Machine& m, double max_size) {
    const std::vector<std::string> names = circuit_names(m);
    std::vector<std::vector<Winding>> out(names.size());
    for (const Coil& c : m.coils) {
        const auto k = static_cast<std::size_t>(
            std::find(names.begin(), names.end(), c.circuit) - names.begin());
        const int nr = std::max(1, static_cast<int>(std::ceil(c.dR / max_size)));
        const int nz = std::max(1, static_cast<int>(std::ceil(c.dZ / max_size)));
        const double turns = c.turns / (nr * nz);
        for (int a = 0; a < nr; ++a) {
            for (int b = 0; b < nz; ++b) {
                out[k].push_back({c.R + c.dR * ((a + 0.5) / nr - 0.5),
                                  c.Z + c.dZ * ((b + 0.5) / nz - 0.5), turns});
            }
        }
    }
    return out;
}

} // namespace shafranov
