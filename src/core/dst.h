#pragma once

#include <span>

namespace shafranov {

// In-place DST-I: X_k = sum_{j=1..n} x_j sin(pi j k / (n + 1)), k = 1..n (0-based storage).
// n + 1 must be a power of two (throws std::invalid_argument otherwise). O(n log n).
// Applying it twice multiplies by (n + 1) / 2, so the inverse is the same transform scaled by 2/(n+1).
void dst1(std::span<double> x);

} // namespace shafranov
