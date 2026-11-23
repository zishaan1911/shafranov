#pragma once

#include "core/equilibrium.h"
#include "core/machine.h"
#include "sim/worker.h"

#include <cstdint>
#include <vector>

namespace shafranov {

struct ViewState {
    Request request;                // what the controls edit
    bool fixed_currents = false;    // "Fixed currents" mode instead of shape control
    std::vector<float> slider;      // kA per circuit (MachineGrid::circuits order), fixed mode
    std::uint64_t seen_id = 0;      // last snapshot the plot caches were built from

    // Plot caches, rebuilt when a new snapshot arrives.
    std::vector<float> heat;        // psi_n clamped to [0, 2], row 0 at the top
    std::vector<double> flux_x, flux_y; // thin contour segments, two points per segment
    std::vector<double> sep_x, sep_y;   // psi_n = 1 segments
};

// Draws the controls and the cross-section into the current ImGui window.
void draw(ViewState& state, Worker& worker, const MachineGrid& mg, const Machine& machine);

} // namespace shafranov
