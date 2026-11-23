#pragma once

#include "core/machine.h"

#include <string>
#include <vector>

namespace shafranov {

// A filament carrying `turns` times its circuit's current.
struct Winding {
    double R, Z, turns;
};

// Circuit names in order of first appearance in m.coils.
std::vector<std::string> circuit_names(const Machine& m);

// Windings of each circuit, in circuit_names order. Each coil's rectangle is split into
// ceil(dR / max_size) x ceil(dZ / max_size) filaments that share its turns equally.
std::vector<std::vector<Winding>> circuit_windings(const Machine& m, double max_size = 0.02);

} // namespace shafranov
