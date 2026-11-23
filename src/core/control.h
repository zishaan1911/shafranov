#pragma once

#include "core/circuits.h"
#include "core/equilibrium.h"

#include <string>
#include <vector>

namespace shafranov {

// Power supply voltage limits (V) per circuit, from Marchioni (EPFL thesis 2024) table 2.2:
// OH1, OH2 1400; E* 651; F* 1250; G 566. Any other circuit (T) gets 0: it is not driven.
std::vector<double> supply_limits(const std::vector<std::string>& circuits);

// Vertical position control on the fast in-vessel G coil: PD on the magnetic axis height,
//   V_G = -(kp (z - z_ref) + kd dz/dt).
// The gains are negative for TCV's G polarity; tuned in the full simulation, they hold a kicked
// kappa 1.4-1.5 plasma to tens of micrometres with tens of volts.
struct VerticalControl {
    double kp = -2e4; // V/m
    double kd = -20.0; // V s/m
    double z_ref = 0.0;
};

// Coil current control for every circuit except G (vertical) and undriven ones:
//   V = R I_ref + M_cc (I_ref - I) / tau
// with M_cc the inductance block of those circuits, so each current relaxes to its reference
// with time constant tau, decoupled from the others.
struct CurrentControl {
    double tau = 5e-3; // s
};

class Controller {
public:
    Controller(const CircuitModel& model, const std::vector<std::string>& circuits);

    // Circuit voltages for the next step, clipped to supply_limits. now: present currents;
    // reference: one target current per circuit; z, z_prev: axis height now and one step earlier.
    // With vertical_on false the G circuit only gets R I (no feedback).
    std::vector<double> voltages(const Currents& now, const std::vector<double>& reference, double z,
                                 double z_prev, double dt) const;

    VerticalControl vertical;
    CurrentControl current;
    bool vertical_on = true;

private:
    const CircuitModel& model_;
    std::vector<double> limits_;
    int g_ = -1;                     // index of the G circuit (-1 if absent)
    std::vector<int> current_controlled; // circuits under current control
};

} // namespace shafranov
