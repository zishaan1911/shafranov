#include "core/poisson.h"

#include "core/constants.h"
#include "core/dst.h"

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace shafranov {

void solve_poisson(const Grid& g, const std::vector<double>& rhs, std::vector<double>& psi) {
    const std::size_t total = static_cast<std::size_t>(g.size());
    if (psi.size() != total || rhs.size() != total) {
        throw std::invalid_argument("solve_poisson: psi or rhs size does not match grid");
    }
    const int m1 = g.nZ - 1;
    if ((m1 & (m1 - 1)) != 0) {
        throw std::invalid_argument("solve_poisson: nZ - 1 must be a power of two");
    }

    const std::size_t n = static_cast<std::size_t>(g.nZ - 2); // interior points per column
    const std::size_t m = static_cast<std::size_t>(g.nR - 2); // interior columns
    const double dR = g.dR();
    const double dZ = g.dZ();
    const double invR2 = 1.0 / (dR * dR);
    const double invZ2 = 1.0 / (dZ * dZ);

    std::vector<double> a(m);
    std::vector<double> c(m);
    for (std::size_t i = 0; i < m; ++i) {
        const double R = g.R(static_cast<int>(i) + 1);
        a[i] = invR2 + 1.0 / (2.0 * R * dR);
        c[i] = invR2 - 1.0 / (2.0 * R * dR);
    }

    // b holds the interior right-hand side, column i (0-based interior) at [i * n, (i + 1) * n).
    std::vector<double> b(m * n);
    for (std::size_t i = 0; i < m; ++i) {
        const int gi = static_cast<int>(i) + 1;
        for (std::size_t j = 0; j < n; ++j) {
            const int gj = static_cast<int>(j) + 1;
            double v = rhs[static_cast<std::size_t>(g.idx(gi, gj))];
            if (i == 0) {
                v -= a[i] * psi[static_cast<std::size_t>(g.idx(0, gj))];
            }
            if (i == m - 1) {
                v -= c[i] * psi[static_cast<std::size_t>(g.idx(g.nR - 1, gj))];
            }
            if (j == 0) {
                v -= invZ2 * psi[static_cast<std::size_t>(g.idx(gi, 0))];
            }
            if (j == n - 1) {
                v -= invZ2 * psi[static_cast<std::size_t>(g.idx(gi, g.nZ - 1))];
            }
            b[i * n + j] = v;
        }
    }

    for (std::size_t i = 0; i < m; ++i) {
        dst1(std::span<double>(b.data() + i * n, n));
    }

    // One tridiagonal solve in R per Z mode (Thomas algorithm).
    std::vector<double> cp(m);
    std::vector<double> dp(m);
    for (std::size_t k = 0; k < n; ++k) {
        const double lambda =
            (2.0 * std::cos(pi * static_cast<double>(k + 1) / static_cast<double>(n + 1)) - 2.0) *
            invZ2;
        const double diag = -2.0 * invR2 + lambda;

        cp[0] = c[0] / diag;
        dp[0] = b[k] / diag;
        for (std::size_t i = 1; i < m; ++i) {
            const double denom = diag - a[i] * cp[i - 1];
            cp[i] = c[i] / denom;
            dp[i] = (b[i * n + k] - a[i] * dp[i - 1]) / denom;
        }
        b[(m - 1) * n + k] = dp[m - 1];
        for (std::size_t i = m - 1; i-- > 0;) {
            b[i * n + k] = dp[i] - cp[i] * b[(i + 1) * n + k];
        }
    }

    const double scale = 2.0 / static_cast<double>(n + 1);
    for (std::size_t i = 0; i < m; ++i) {
        std::span<double> col(b.data() + i * n, n);
        dst1(col);
        for (std::size_t j = 0; j < n; ++j) {
            psi[static_cast<std::size_t>(g.idx(static_cast<int>(i) + 1, static_cast<int>(j) + 1))] =
                scale * col[j];
        }
    }
}

} // namespace shafranov
