#pragma once

namespace shafranov {

struct BField {
    double R;
    double Z;
};

// Flux per radian at (R, Z) from a 1 A filament at (Rc, Zc). Symmetric in the two points.
double greens_psi(double Rc, double Zc, double R, double Z);

// Field at (R, Z) from a 1 A filament at (Rc, Zc). Requires R > 0.
BField greens_b(double Rc, double Zc, double R, double Z);

// Mutual inductance between two filaments (H): 2 pi psi.
double mutual_inductance(double R1, double Z1, double R2, double Z2);

} // namespace shafranov
