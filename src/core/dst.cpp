#include "core/dst.h"

#include "core/constants.h"

#include <cstddef>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

namespace shafranov {

namespace {

using Complex = std::complex<double>;

// Iterative radix-2 Cooley-Tukey FFT (forward, e^{-i...}); a.size() must be a power of two.
void fft(std::vector<Complex>& a) {
    const std::size_t n = a.size();

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len / 2;
        for (std::size_t start = 0; start < n; start += len) {
            for (std::size_t k = 0; k < half; ++k) {
                // Direct evaluation per index keeps twiddle error from accumulating.
                const Complex w = std::polar(1.0, -2.0 * pi * static_cast<double>(k) /
                                                      static_cast<double>(len));
                const Complex u = a[start + k];
                const Complex v = a[start + k + half] * w;
                a[start + k] = u + v;
                a[start + k + half] = u - v;
            }
        }
    }
}

} // namespace

void dst1(std::span<double> x) {
    const std::size_t n = x.size();
    const std::size_t m = n + 1;
    if ((m & (m - 1)) != 0) {
        throw std::invalid_argument("dst1: n + 1 must be a power of two");
    }
    if (n == 0) {
        return;
    }

    // Odd extension y = [0, x1..xn, 0, -xn..-x1].
    std::vector<Complex> y(2 * m, Complex(0.0, 0.0));
    for (std::size_t j = 1; j <= n; ++j) {
        y[j] = x[j - 1];
        y[2 * m - j] = -x[j - 1];
    }

    fft(y);

    for (std::size_t k = 1; k <= n; ++k) {
        x[k - 1] = -y[k].imag() / 2.0;
    }
}

} // namespace shafranov
