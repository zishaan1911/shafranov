#pragma once

#include "core/contour.h"
#include "core/equilibrium.h"
#include "core/inverse.h"

#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace shafranov {

// What the user asks for: a plasma shape (the inverse solve picks the E and F coil currents), or,
// with fixed_currents set, given circuit currents (a forward solve: no shape control).
struct Request {
    double Ip = 2e5;      // A
    double beta0 = 0.3;
    double R0 = 0.88;     // boundary centre (m)
    double a = 0.22;      // minor radius (m)
    double kappa = 1.4;
    double delta = 0.3;
    bool diverted = true; // lower single null with an X-point target
    double xR = 0.814, xZ = -0.339;
    std::optional<std::vector<double>> fixed_currents; // per circuit, MachineGrid::circuits order
};

// Isoflux targets for a request: 12 points on the Miller boundary
//   R = R0 + a cos(t + delta sin t), Z = kappa a sin t,  t = k 2 pi / 12,
// dropping the bottom point (k = 9) when diverted, paired with the X-point (diverted) or with the
// outboard midplane point t = 0 (limited). Fixed circuits: OH1, OH2, G, T.
ShapeTargets targets_for(const Request& r);

struct Snapshot {
    std::uint64_t id = 0;   // increases with every finished solve
    Request request;
    Equilibrium eq;
    Currents currents;
    std::optional<ShapeMetrics> metrics;
    std::string error;      // empty on success; the solver's message if it failed (e.g. lost plasma)
    double solve_ms = 0.0;
};

// Solves on a background thread. request() replaces any pending request (only the newest matters,
// so dragging a slider never queues work). Each solve warm-starts from the last successful one.
class Worker {
public:
    explicit Worker(const MachineGrid& mg);
    ~Worker();
    Worker(const Worker&) = delete;
    Worker& operator=(const Worker&) = delete;

    void request(const Request& r);
    std::shared_ptr<const Snapshot> latest() const; // null until the first solve finishes
    bool busy() const;                              // a solve is running or pending

private:
    void run();

    const MachineGrid& mg_;
    mutable std::mutex mutex_;
    std::condition_variable wake_;
    std::optional<Request> pending_;
    bool solving_ = false;
    bool stop_ = false;
    std::shared_ptr<const Snapshot> latest_;
    std::thread thread_; // last member: starts after everything above is constructed
};

} // namespace shafranov
