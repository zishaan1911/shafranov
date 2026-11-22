#include "core/profiles.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "core/constants.h"

namespace shafranov {

namespace {

void validate(const Profile& p) {
    if (!std::isfinite(p.Ip)) {
        throw std::invalid_argument("Profile: Ip must be finite");
    }
    if (!(p.beta0 >= 0.0 && p.beta0 <= 1.0)) {
        throw std::invalid_argument("Profile: beta0 must be in [0, 1]");
    }
    if (!(p.alpha_m > 0.0) || !(p.alpha_n > 0.0)) {
        throw std::invalid_argument("Profile: alpha_m and alpha_n must be positive");
    }
    if (!(p.R0 > 0.0)) {
        throw std::invalid_argument("Profile: R0 must be positive");
    }
}

double shape(const Profile& p, double psi_n) {
    const double x = std::clamp(psi_n, 0.0, 1.0);
    return std::pow(1.0 - std::pow(x, p.alpha_m), p.alpha_n);
}

} // namespace

double plasma_current(const Grid& g, const std::vector<double>& psi, double psi_axis,
                      double psi_bndry, const std::vector<char>& inside, const Profile& p,
                      std::vector<double>& jphi) {
    const auto n = static_cast<std::size_t>(g.size());
    if (psi.size() != n || inside.size() != n) {
        throw std::invalid_argument("plasma_current: psi and inside sizes must match grid");
    }
    if (psi_bndry == psi_axis) {
        throw std::invalid_argument("plasma_current: psi_bndry must differ from psi_axis");
    }
    validate(p);

    jphi.assign(n, 0.0);
    double sum = 0.0;
    for (int i = 0; i < g.nR; ++i) {
        const double R = g.R(i);
        const double radial = p.beta0 * R / p.R0 + (1.0 - p.beta0) * p.R0 / R;
        for (int j = 0; j < g.nZ; ++j) {
            const auto k = static_cast<std::size_t>(g.idx(i, j));
            if (inside[k] == 0) {
                continue;
            }
            const double psi_n = (psi[k] - psi_axis) / (psi_bndry - psi_axis);
            jphi[k] = radial * shape(p, psi_n);
            sum += jphi[k];
        }
    }

    const double total = sum * g.dR() * g.dZ();
    if (total == 0.0) {
        throw std::runtime_error("plasma_current: unscaled current integral is zero");
    }
    const double lambda = p.Ip / total;
    for (double& v : jphi) {
        v *= lambda;
    }
    return lambda;
}

double pprime(const Profile& p, double lambda, double psi_n) {
    return lambda * p.beta0 / p.R0 * shape(p, psi_n);
}

double ffprime(const Profile& p, double lambda, double psi_n) {
    return mu0 * lambda * (1.0 - p.beta0) * p.R0 * shape(p, psi_n);
}

} // namespace shafranov
