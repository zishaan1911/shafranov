#include "core/control.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <stdexcept>

namespace shafranov {

std::vector<double> supply_limits(const std::vector<std::string>& circuits) {
    std::vector<double> limits;
    limits.reserve(circuits.size());
    for (const std::string& name : circuits) {
        double limit = 0.0;
        if (name == "OH1" || name == "OH2") {
            limit = 1400.0;
        } else if (name == "G") {
            limit = 566.0;
        } else if (name.size() >= 2 && std::isdigit(static_cast<unsigned char>(name[1])) != 0) {
            if (name[0] == 'E') {
                limit = 651.0;
            } else if (name[0] == 'F') {
                limit = 1250.0;
            }
        }
        limits.push_back(limit);
    }
    return limits;
}

Controller::Controller(const CircuitModel& model, const std::vector<std::string>& circuits)
    : model_(model), limits_(supply_limits(circuits)) {
    for (std::size_t c = 0; c < circuits.size(); ++c) {
        if (circuits[c] == "G") {
            g_ = static_cast<int>(c);
        } else if (limits_[c] > 0.0) {
            current_controlled.push_back(static_cast<int>(c));
        }
    }
}

std::vector<double> Controller::voltages(const Currents& now, const std::vector<double>& reference,
                                         double z, double z_prev, double dt) const {
    if (now.circuits.size() != limits_.size() || reference.size() != limits_.size()) {
        throw std::invalid_argument("Controller::voltages: one entry per circuit required");
    }
    std::vector<double> v(limits_.size(), 0.0);

    for (const int c : current_controlled) {
        const auto uc = static_cast<std::size_t>(c);
        double volts = model_.R[c] * reference[uc];
        for (const int d : current_controlled) {
            const auto ud = static_cast<std::size_t>(d);
            volts += model_.M(c, d) * (reference[ud] - now.circuits[ud]) / current.tau;
        }
        v[uc] = volts;
    }

    if (g_ >= 0) {
        const auto ug = static_cast<std::size_t>(g_);
        double volts = model_.R[g_] * now.circuits[ug];
        if (vertical_on) {
            volts -= vertical.kp * (z - vertical.z_ref) + vertical.kd * (z - z_prev) / dt;
        }
        v[ug] = volts;
    }

    for (std::size_t c = 0; c < v.size(); ++c) {
        v[c] = std::clamp(v[c], -limits_[c], limits_[c]);
    }
    return v;
}

} // namespace shafranov
