#pragma once

#include "core/control.h"
#include "core/dynamics.h"
#include "sim/worker.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace shafranov {

// A scripted discharge: plasma shapes (and Ip) at given times. Between waypoints the coil current
// references and Ip are interpolated linearly; after the last one they hold.
struct Waypoint {
    double t; // s
    Request shape;
};

// TCV-like demo: Ip 120 kA at kappa 1.25 rising to 200 kA at kappa 1.45 by 20 ms, held to 40 ms,
// lower single null throughout.
std::vector<Waypoint> demo_scenario();

// Reference circuit currents for each waypoint, from the inverse solve (each warm-started from the
// previous). Also returns the starting equilibrium (waypoint 0). Throws std::runtime_error if a
// waypoint does not converge.
struct Plan {
    std::vector<Waypoint> waypoints;
    std::vector<std::vector<double>> currents; // per waypoint, per circuit
    std::vector<double> axis_z;                // planned magnetic axis height per waypoint (m)
    Equilibrium start;
};
Plan make_plan(const MachineGrid& mg, const std::vector<Waypoint>& waypoints);

// Interpolated references at time t.
std::vector<double> reference_currents(const Plan& plan, double t);
double reference_ip(const Plan& plan, double t);
// Planned axis height at t: the vertical controller's target, so it does not fight the coils.
double reference_z(const Plan& plan, double t);

// One recorded instant of the discharge.
struct Sample {
    double t, R_axis, Z_axis, Ip, I_G, V_G;
    bool vertical_on;
};

enum class DischargeStatus { Paused, Running, Finished, Lost };

struct DischargeView {
    DischargeStatus status = DischargeStatus::Paused;
    double t = 0.0;
    Equilibrium eq;
    Currents currents;
    std::vector<Sample> history;
    std::string message; // why the plasma was lost, if it was
};

// Runs the discharge on a background thread with a Controller, one Simulation step (dt = 1e-4 s)
// at a time, so the window can watch it. The plan is computed in the constructor.
class Discharge {
public:
    Discharge(const MachineGrid& mg, const CircuitModel& model, std::vector<Waypoint> waypoints,
              double end_time = 0.04);
    ~Discharge();
    Discharge(const Discharge&) = delete;
    Discharge& operator=(const Discharge&) = delete;

    void run();                      // start or resume
    void pause();
    void reset();                    // back to t = 0, paused
    void set_vertical_control(bool on);
    void kick(double volts = 20.0);  // one step of extra G voltage: a vertical nudge
    std::shared_ptr<const DischargeView> view() const;

private:
    void loop();

    const MachineGrid& mg_;
    const CircuitModel& model_;
    Plan plan_;
    double end_time_;
    mutable std::mutex mutex_;
    bool running_ = false, stop_ = false, reset_ = false, vertical_ = true;
    double kick_ = 0.0;
    std::shared_ptr<const DischargeView> view_;
    std::thread thread_; // last: starts after everything above
};

} // namespace shafranov
