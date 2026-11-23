#pragma once

#include "core/equilibrium.h"
#include "core/machine.h"
#include "sim/discharge.h"
#include "sim/worker.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace shafranov {

// Cross-section plot data, rebuilt when a new equilibrium arrives.
struct PlotCache {
    std::vector<float> heat;        // psi_n clamped to [0, 2], row 0 at the top
    std::vector<double> flux_x, flux_y; // thin contour segments, two points per segment
    std::vector<double> sep_x, sep_y;   // psi_n = 1 segments
};

struct ViewState {
    Request request;                // what the controls edit
    bool fixed_currents = false;    // "Fixed currents" mode instead of shape control
    std::vector<float> slider;      // kA per circuit (MachineGrid::circuits order), fixed mode
    std::uint64_t seen_id = 0;      // last snapshot the plot caches were built from
    PlotCache eq_cache;             // Equilibrium mode, from the worker snapshot

    bool discharge_mode = false;    // "Discharge" mode instead of "Equilibrium"
    bool vertical_control = true;   // mirrors Discharge::set_vertical_control
    std::shared_ptr<const DischargeView> seen_view; // view dis_cache was built from
    PlotCache dis_cache;            // Discharge mode, from the discharge view
    double t_lo = 0.0, t_hi = 40.0; // linked x limits of the time traces (ms)
};

// Draws the controls and the cross-section into the current ImGui window.
void draw(ViewState& state, Worker& worker, Discharge& discharge, const MachineGrid& mg,
          const Machine& machine);

} // namespace shafranov
