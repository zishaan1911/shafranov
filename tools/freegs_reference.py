"""Generate the FreeGS reference equilibrium that tests/test_freegs_reference.cpp compares against.

Both codes solve the same problem: TCV coils from data/tcv/coils.csv, split into the same 2 cm filaments
as src/core/coils.cpp, no vessel current, the same grid, the same fixed-beta0 current profile and the same
shape targets.

Usage: python tools/freegs_reference.py
Needs FreeGS newer than 0.8.2 for NumPy 2 (tested with git commit 5b41fe556):
  pip install "git+https://github.com/freegs-plasma/freegs@5b41fe556"
Writes tests/data/freegs_tcv_lsn.geqdsk and tests/data/freegs_tcv_lsn.csv (circuit currents, X-point).
"""

import math
import pathlib

import freegs
import numpy as np
from freegs import critical, geqdsk
from freegs.machine import Circuit, Coil, Machine

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "tests" / "data"

# Shared with tests/test_inverse.cpp and tests/test_freegs_reference.cpp.
GRID = dict(Rmin=0.55, Rmax=1.20, Zmin=-0.80, Zmax=0.80, nx=65, ny=129)
IP, BETA0, R0 = 2e5, 0.3, 0.88
XPOINT = (0.814, -0.339)
FIXED = {"OH1", "OH2", "G", "T"}
GAMMA = 1e-12  # Tikhonov weight on currents (A^-2), as ShapeTargets::gamma


def boundary_points():
    pts = []
    for k in range(12):
        if k == 9:  # the bottom, replaced by the X-point
            continue
        t = k * 2 * math.pi / 12
        pts.append((0.88 + 0.22 * math.cos(t + 0.3 * math.sin(t)), 1.4 * 0.22 * math.sin(t)))
    return pts


def read_coils():
    rows = [line for line in (ROOT / "data" / "tcv" / "coils.csv").read_text().splitlines()
            if line and not line.startswith("#")]
    head = rows[0].split(",")
    return [dict(zip(head, r.split(","))) for r in rows[1:]]


def machine():
    """Same split as circuit_windings(m, 0.02)."""
    circuits, order = {}, []
    for c in read_coils():
        R, Z, dR, dZ, turns = (float(c[k]) for k in ("R", "Z", "dR", "dZ", "turns"))
        nr, nz = max(1, math.ceil(dR / 0.02)), max(1, math.ceil(dZ / 0.02))
        name = c["circuit"]
        if name not in circuits:
            circuits[name] = []
            order.append(name)
        for a in range(nr):
            for b in range(nz):
                circuits[name].append((f"{c['name']}_{a}_{b}",
                                       Coil(R + dR * ((a + 0.5) / nr - 0.5),
                                            Z + dZ * ((b + 0.5) / nz - 0.5)),
                                       turns / (nr * nz)))
    return Machine([(n, Circuit(circuits[n], current=0.0, control=n not in FIXED)) for n in order]), order


class FixedBeta0(freegs.jtor.ConstrainPaxisIp):
    """J = L (beta0 R/R0 + (1-beta0) R0/R) shape, with beta0 fixed and L set by Ip (as in profiles.cpp)."""

    def Jtor(self, R, Z, psi, psi_bndry=None):
        opt, xpt = critical.find_critical(R, Z, psi)
        if not opt:
            raise ValueError("No O-points found!")
        psi_axis = opt[0][2]
        if psi_bndry is not None:
            mask = critical.core_mask(R, Z, psi, opt, xpt, psi_bndry)
        elif xpt:
            psi_bndry = xpt[0][2]
            mask = critical.core_mask(R, Z, psi, opt, xpt)
        else:
            psi_bndry, mask = psi[0, 0], None
        dR, dZ = R[1, 0] - R[0, 0], Z[0, 1] - Z[0, 0]
        psi_norm = (psi - psi_axis) / (psi_bndry - psi_axis)
        shape = (1.0 - np.clip(psi_norm, 0.0, 1.0) ** self.alpha_m) ** self.alpha_n
        if mask is not None:
            shape *= mask
        radial = BETA0 * R / R0 + (1 - BETA0) * R0 / R
        self.Beta0 = BETA0
        self.L = IP / (np.sum(radial * shape) * dR * dZ)
        self.psi_axis, self.psi_bndry = psi_axis, psi_bndry
        return self.L * radial * shape


class AbsoluteTikhonov:
    """Like freegs.control.constrain, but regularises the absolute currents, not the change, with weight
    GAMMA (no square): x = argmin |A x - b|^2 + GAMMA |x|^2. This is what src/core/inverse.cpp does."""

    def __init__(self, xpoints, isoflux):
        self.xpoints, self.isoflux = xpoints, isoflux

    def __call__(self, eq):
        tokamak = eq.getMachine()
        A, b = [], []
        for R, Z in self.xpoints:
            A.append(tokamak.controlBr(R, Z))
            b.append(-eq.Br(R, Z))
            A.append(tokamak.controlBz(R, Z))
            b.append(-eq.Bz(R, Z))
        for r1, z1, r2, z2 in self.isoflux:
            A.append([c1 - c2 for c1, c2 in zip(tokamak.controlPsi(r1, z1), tokamak.controlPsi(r2, z2))])
            b.append(eq.psiRZ(r2, z2) - eq.psiRZ(r1, z1))
        A, b = np.array(A), np.reshape(np.array(b, dtype=float), (-1,))
        x0 = np.array(tokamak.controlCurrents())
        x = np.linalg.solve(A.T @ A + GAMMA * np.eye(len(x0)), A.T @ (b + A @ x0))
        tokamak.controlAdjust(x - x0)


def main():
    tokamak, order = machine()
    eq = freegs.Equilibrium(tokamak=tokamak, boundary=freegs.boundary.freeBoundary, **GRID)
    profiles = FixedBeta0(eq, paxis=1e3, Ip=IP, fvac=0.88 * 1.4, Raxis=R0)
    pts = boundary_points()
    constrain = AbsoluteTikhonov([XPOINT], [(XPOINT[0], XPOINT[1], R, Z) for R, Z in pts])
    freegs.solve(eq, profiles, constrain, rtol=1e-6, maxits=200)

    _, xpt = critical.find_critical(eq.R, eq.Z, eq.psi())
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "freegs_tcv_lsn.geqdsk", "w") as f:
        geqdsk.write(eq, f)
    rows = ["quantity,value"] + [f"{n},{float(tokamak[n].current)!r}" for n in order]
    rows += [f"xpoint_R,{float(xpt[0][0])!r}", f"xpoint_Z,{float(xpt[0][1])!r}"]
    (OUT / "freegs_tcv_lsn.csv").write_text("\n".join(rows) + "\n", newline="\n")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
