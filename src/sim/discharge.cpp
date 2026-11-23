#include "sim/discharge.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace shafranov {

namespace {

constexpr double step_dt = 1e-4;

Profile profile_of(const Request& r, double Ip) {
    return Profile{Ip, r.beta0, 1.0, 2.0, r.R0};
}

// Position in the plan: the waypoint index below t and the fraction towards the next one.
struct Bracket {
    std::size_t lo, hi;
    double w;
};

Bracket locate(const Plan& plan, double t) {
    const auto& wp = plan.waypoints;
    if (t <= wp.front().t) {
        return {0, 0, 0.0};
    }
    if (t >= wp.back().t) {
        return {wp.size() - 1, wp.size() - 1, 0.0};
    }
    std::size_t hi = 1;
    while (wp[hi].t < t) {
        ++hi;
    }
    return {hi - 1, hi, (t - wp[hi - 1].t) / (wp[hi].t - wp[hi - 1].t)};
}

} // namespace

std::vector<Waypoint> demo_scenario() {
    Request start;
    start.Ip = 120e3;
    start.kappa = 1.25;
    Request flat;
    flat.Ip = 200e3;
    flat.kappa = 1.45;
    return {{0.0, start}, {0.02, flat}, {0.04, flat}};
}

Plan make_plan(const MachineGrid& mg, const std::vector<Waypoint>& waypoints) {
    if (waypoints.empty()) {
        throw std::invalid_argument("a discharge needs at least one waypoint");
    }
    for (std::size_t k = 1; k < waypoints.size(); ++k) {
        if (!(waypoints[k].t > waypoints[k - 1].t)) {
            throw std::invalid_argument("waypoint times must be strictly increasing");
        }
    }

    Plan plan;
    plan.waypoints = waypoints;
    Currents currents;
    currents.circuits.assign(mg.circuits.size(), 0.0);
    std::vector<double> jphi;
    for (std::size_t k = 0; k < waypoints.size(); ++k) {
        const Request& r = waypoints[k].shape;
        InitialPlasma initial{r.R0, 0.0, 0.15};
        initial.jphi = jphi;
        Equilibrium eq = solve_equilibrium(mg, currents, profile_of(r, r.Ip), initial,
                                           SolveOptions{}, shape_feedback(mg, targets_for(r)));
        if (!eq.converged) {
            throw std::runtime_error("waypoint " + std::to_string(k) + " (t = " +
                                     std::to_string(waypoints[k].t) + " s) did not converge");
        }
        jphi = eq.jphi;
        plan.currents.push_back(currents.circuits);
        plan.axis_z.push_back(eq.topology.axis.Z);
        if (k == 0) {
            plan.start = std::move(eq);
        }
    }
    return plan;
}

std::vector<double> reference_currents(const Plan& plan, double t) {
    const Bracket s = locate(plan, t);
    std::vector<double> out = plan.currents[s.lo];
    for (std::size_t c = 0; c < out.size(); ++c) {
        out[c] += s.w * (plan.currents[s.hi][c] - plan.currents[s.lo][c]);
    }
    return out;
}

double reference_ip(const Plan& plan, double t) {
    const Bracket s = locate(plan, t);
    const double a = plan.waypoints[s.lo].shape.Ip;
    const double b = plan.waypoints[s.hi].shape.Ip;
    return a + s.w * (b - a);
}

double reference_z(const Plan& plan, double t) {
    const Bracket s = locate(plan, t);
    return plan.axis_z[s.lo] + s.w * (plan.axis_z[s.hi] - plan.axis_z[s.lo]);
}

Discharge::Discharge(const MachineGrid& mg, const CircuitModel& model,
                     std::vector<Waypoint> waypoints, double end_time)
    : mg_(mg), model_(model), plan_(make_plan(mg, waypoints)), end_time_(end_time) {
    auto first = std::make_shared<DischargeView>();
    first->eq = plan_.start;
    first->currents.circuits = plan_.currents.front();
    first->currents.vessel.assign(static_cast<std::size_t>(model_.n_vessel), 0.0);
    view_ = std::move(first);
    thread_ = std::thread([this] { loop(); });
}

Discharge::~Discharge() {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
    }
    thread_.join();
}

void Discharge::run() {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        running_ = true;
    }
}

void Discharge::pause() {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        running_ = false;
    }
}

void Discharge::reset() {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        reset_ = true;
    }
}

void Discharge::set_vertical_control(bool on) {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        vertical_ = on;
    }
}

void Discharge::kick(double volts) {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        kick_ = volts;
    }
}

std::shared_ptr<const DischargeView> Discharge::view() const {
    const std::lock_guard<std::mutex> lock(mutex_);
    return view_;
}

void Discharge::loop() {
    const Request& base = plan_.waypoints.front().shape;
    const auto g_index = static_cast<std::size_t>(
        std::find(mg_.circuits.begin(), mg_.circuits.end(), "G") - mg_.circuits.begin());

    std::unique_ptr<Simulation> sim;
    std::unique_ptr<Controller> controller;
    std::vector<Sample> history;
    double z_prev = 0.0;
    bool over = false; // finished or lost: nothing more to step until a reset

    const auto rebuild = [&] {
        Currents start;
        start.circuits = plan_.currents.front();
        start.vessel.assign(static_cast<std::size_t>(model_.n_vessel), 0.0);
        sim = std::make_unique<Simulation>(mg_, model_, profile_of(base, base.Ip), plan_.start,
                                           start, step_dt);
        controller = std::make_unique<Controller>(model_, mg_.circuits);
        // Hold the plasma where its planned equilibrium puts it, not at Z = 0.
        controller->vertical.z_ref = plan_.start.topology.axis.Z;
        history.clear();
        z_prev = plan_.start.topology.axis.Z;
        over = false;
    };
    const auto publish = [&](DischargeStatus status, const std::string& message) {
        auto v = std::make_shared<DischargeView>();
        v->status = status;
        v->t = sim->time();
        v->eq = sim->equilibrium();
        v->currents = sim->currents();
        v->history = history;
        v->message = message;
        const std::lock_guard<std::mutex> lock(mutex_);
        view_ = std::move(v);
    };

    rebuild();

    for (;;) {
        bool vertical = true;
        double kick = 0.0;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            if (stop_) {
                return;
            }
            if (!running_ && !reset_) {
                // discharge.h has no condition variable, so idle by polling.
                lock.unlock();
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                continue;
            }
            if (reset_) {
                lock.unlock();
                rebuild();
                publish(DischargeStatus::Paused, "");
                lock.lock();
                reset_ = false;
                running_ = false;
                continue;
            }
            if (over) {
                running_ = false;
                continue;
            }
            vertical = vertical_;
            kick = kick_;
            kick_ = 0.0;
        }

        const double t = sim->time();
        const double z = sim->equilibrium().topology.axis.Z;
        controller->vertical_on = vertical;
        const std::vector<double> ref = reference_currents(plan_, t);
        controller->vertical.z_ref = reference_z(plan_, t);
        std::vector<double> V =
            controller->voltages(sim->currents(), ref, z, z_prev, sim->dt());
        V[g_index] += kick;
        const double ip = reference_ip(plan_, t + sim->dt());
        sim->set_profile(profile_of(base, ip));

        try {
            sim->step(V);
        } catch (const PlasmaLost& e) {
            over = true;
            publish(DischargeStatus::Lost, e.what());
            const std::lock_guard<std::mutex> lock(mutex_);
            running_ = false;
            continue;
        }

        z_prev = z;
        const CriticalPoint& axis = sim->equilibrium().topology.axis;
        history.push_back({sim->time(), axis.R, axis.Z, ip, sim->currents().circuits[g_index],
                           V[g_index], vertical});
        if (sim->time() >= end_time_) {
            over = true;
            publish(DischargeStatus::Finished, "");
            const std::lock_guard<std::mutex> lock(mutex_);
            running_ = false;
        } else {
            publish(DischargeStatus::Running, "");
        }
    }
}

} // namespace shafranov
