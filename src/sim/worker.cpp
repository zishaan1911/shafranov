#include "sim/worker.h"

#include "core/constants.h"

#include <chrono>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <utility>

namespace shafranov {

ShapeTargets targets_for(const Request& r) {
    ShapeTargets targets;
    const Point anchor = r.diverted ? Point{r.xR, r.xZ} : Point{r.R0 + r.a, 0.0};
    for (int k = 0; k < 12; ++k) {
        if (r.diverted && k == 9) {
            continue; // the bottom, replaced by the X-point
        }
        if (!r.diverted && k == 0) {
            continue; // the anchor itself
        }
        const double t = k * 2.0 * pi / 12.0;
        const double R = r.R0 + r.a * std::cos(t + r.delta * std::sin(t));
        targets.isoflux.push_back({anchor.R, anchor.Z, R, r.kappa * r.a * std::sin(t)});
    }
    if (r.diverted) {
        targets.xpoints.push_back(anchor);
    }
    targets.fixed = {"OH1", "OH2", "G", "T"};
    return targets;
}

Worker::Worker(const MachineGrid& mg) : mg_(mg), thread_([this] { run(); }) {}

Worker::~Worker() {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
    }
    wake_.notify_all();
    thread_.join();
}

void Worker::request(const Request& r) {
    {
        const std::lock_guard<std::mutex> lock(mutex_);
        pending_ = r;
    }
    wake_.notify_all();
}

std::shared_ptr<const Snapshot> Worker::latest() const {
    const std::lock_guard<std::mutex> lock(mutex_);
    return latest_;
}

bool Worker::busy() const {
    const std::lock_guard<std::mutex> lock(mutex_);
    return pending_.has_value() || solving_;
}

void Worker::run() {
    std::shared_ptr<const Snapshot> last_good; // only touched by this thread
    std::uint64_t id = 0;

    for (;;) {
        Request r;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            wake_.wait(lock, [this] { return pending_.has_value() || stop_; });
            if (stop_) {
                return;
            }
            r = *pending_;
            pending_.reset();
            solving_ = true;
        }

        Snapshot snap;
        snap.request = r;
        const auto start = std::chrono::steady_clock::now();
        bool ok = false;
        try {
            Currents currents;
            if (r.fixed_currents) {
                if (r.fixed_currents->size() != mg_.circuits.size()) {
                    throw std::invalid_argument("fixed_currents needs one value per circuit");
                }
                currents.circuits = *r.fixed_currents;
            } else if (last_good) {
                currents.circuits = last_good->currents.circuits;
            } else {
                currents.circuits.assign(mg_.circuits.size(), 0.0);
            }

            InitialPlasma initial{r.R0, 0.0, 0.7 * r.a};
            if (last_good) {
                initial.jphi = last_good->eq.jphi;
            }
            const Profile profile{r.Ip, r.beta0, 1.0, 2.0, r.R0};
            Feedback feedback;
            if (!r.fixed_currents) {
                feedback = shape_feedback(mg_, targets_for(r));
            }

            snap.eq = solve_equilibrium(mg_, currents, profile, initial, SolveOptions{}, feedback);
            snap.currents = std::move(currents);
            ok = true;
            try {
                snap.metrics = shape_metrics(mg_.grid, snap.eq.psi, snap.eq.topology);
            } catch (const std::exception&) {
                snap.metrics.reset();
            }
            if (!snap.eq.converged) {
                snap.error = "did not converge in " + std::to_string(snap.eq.iterations) +
                             " iterations";
            }
        } catch (const std::exception& e) {
            snap.error = e.what();
            snap.metrics.reset();
            if (last_good) {
                snap.eq = last_good->eq;
                snap.currents = last_good->currents;
            }
        }
        snap.solve_ms = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - start)
                            .count();
        snap.id = ++id;

        auto published = std::make_shared<const Snapshot>(std::move(snap));
        if (ok && published->eq.converged) { // an unconverged state is a poor warm start
            last_good = published;
        }
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            latest_ = std::move(published);
            solving_ = false;
        }
    }
}

} // namespace shafranov
