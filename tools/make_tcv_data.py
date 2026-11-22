"""Regenerate data/tcv/*.csv from the published sources.

Sources (pinned to exact commits):
  RZIP TCV structure model (J Wainwright, A Sharma; GPL-2 code, data used with attribution):
    https://github.com/atisharma/rzip  rzip/tcv-tiles/tcv_structure.mat
    Same quantities as the TCV MDSplus STATIC tree (R_V, Z_V, RES_V, coil cross-sections).
  Coil table and circuits: arXiv:1301.7617, appendix tables 2 and 3.
  First wall: LIUQE equilibrium of TCV shot 65402 at t = 1 s, from the TCV-X21 dataset
    (Oliveira, Body et al., Nucl. Fusion 62 (2022) 096001; CC-BY-4.0):
    https://github.com/SPCData/TCV-X21  1.experimental_data/reference_scenario/65402_t1.eqdsk

Usage: python tools/make_tcv_data.py   (needs numpy and scipy)
"""

import io
import math
import pathlib
import urllib.request

import numpy as np
import scipy.io

RZIP = ("https://raw.githubusercontent.com/atisharma/rzip/"
        "014dce1dd58e61e4178dbae07db6c5580ec79506/rzip/tcv-tiles/tcv_structure.mat")
X21 = ("https://raw.githubusercontent.com/SPCData/TCV-X21/"
       "9a0a5eb13e2f0ff34494c31a66a3d32b042ddd7d/1.experimental_data/reference_scenario/"
       "65402_t1.eqdsk")
OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "tcv"


def fetch(url):
    with urllib.request.urlopen(url) as r:
        return r.read()


def num(x):
    # The .mat stores float32; print the shortest decimal that round-trips it.
    return str(np.float32(x))


def coils(m):
    names = (["A1", "B1", "B2", "C1", "C2", "D1", "D2"]
             + [f"E{k}" for k in range(1, 9)] + [f"F{k}" for k in range(1, 9)]
             + [f"G{k}" for k in range(1, 7)] + ["T1", "T2", "T3"])
    circuits = (["OH1"] + ["OH2"] * 6 + [f"E{k}" for k in range(1, 9)]
                + [f"F{k}" for k in range(1, 9)] + ["G"] * 6 + ["T"] * 3)
    # arXiv:1301.7617 table 3: lower G windings +1, upper -1; T busbar weights.
    weights = [1.0] * 23 + [1, 1, 1, -1, -1, -1] + [-26 / 68, -42 / 68, 1.0]
    rho_cu = float(m["fil_re"])
    rows = ["name,circuit,R,Z,dR,dZ,turns,resistance"]
    for k, name in enumerate(names):
        R, Z = float(m["coil_r"][k]), float(m["coil_z"][k])
        w, h = float(m["coil_w"][k]), float(m["coil_h"][k])
        if math.isnan(w):  # G windings: the .mat leaves their size blank; table 2 gives 5 mm
            w = h = 0.005
        n = int(m["nt_c"][k])
        turns = n * weights[k]
        # Series resistance of n turns of conductor cross-section cs around the torus.
        res = rho_cu * n * 2 * math.pi * R / float(m["coil_cs"][k])
        rows.append(f"{name},{circuits[k]},{num(R)},{num(Z)},{num(w)},{num(h)},"
                    f"{turns:.6g},{res:.6g}")
    return rows


def vessel(m):
    rows = ["R,Z,dR,dZ,resistance"]
    for R, Z, w, h, res in zip(m["ves_r"], m["ves_z"], m["ves_w"], m["ves_h"], m["res_v"]):
        rows.append(f"{num(R)},{num(Z)},{num(w)},{num(h)},{res:.6g}")
    return rows


def limiter(text):
    lines = text.splitlines()
    nw, nh = map(int, lines[0].split()[-2:])
    pos = 1

    def take(n):
        nonlocal pos
        out = []
        while len(out) < n:
            line = lines[pos].rstrip()
            pos += 1
            out += [float(line[k:k + 16]) for k in range(0, len(line), 16)]
        return out

    take(20)
    for n in [nw, nw, nw, nw, nw * nh, nw]:
        take(n)
    nbbbs, limitr = int(lines[pos][:5]), int(lines[pos][5:10])
    pos += 1
    take(2 * nbbbs)
    lim = take(2 * limitr)
    return ["R,Z"] + [f"{r:.6g},{z:.6g}" for r, z in zip(lim[0::2], lim[1::2])]


HEADERS = {
    "coils.csv": [
        "# TCV poloidal field coils.",
        "# Geometry and turns: RZIP tcv_structure.mat (coil_r/z/w/h, nt_c), identical to",
        "#   arXiv:1301.7617 table 2. G winding size from table 2.",
        "# Circuits and signed turns: arXiv:1301.7617 table 3.",
        "# Resistance: rho_Cu * |turns| * 2 pi R / conductor cross-section, with rho_Cu = fil_re and",
        "#   cross-section = coil_cs from the same .mat. Matches Marchioni (EPFL thesis 2024) table",
        "#   2.2: OH1 15 mOhm, F 39.1 mOhm, E 11.8 mOhm (here 12.5).",
    ],
    "vessel.csv": [
        "# TCV vacuum vessel: 256 toroidal filaments, stainless steel (rho = 7.5e-7 Ohm m).",
        "# RZIP tcv_structure.mat (ves_r, ves_z, ves_w, ves_h, res_v), the MDSplus STATIC",
        "#   R_V, Z_V, RES_V model. Parallel toroidal resistance 55.7 uOhm (Hofmann 1996: 49 uOhm).",
    ],
    "limiter.csv": [
        "# TCV first wall, as used by LIUQE for shot 65402 at t = 1 s.",
        "# TCV-X21 dataset, 65402_t1.eqdsk (CC-BY-4.0): Oliveira, Body et al.,",
        "#   Nucl. Fusion 62 (2022) 096001.",
    ],
}


def main():
    m = scipy.io.loadmat(io.BytesIO(fetch(RZIP)), squeeze_me=True)
    eqdsk = fetch(X21).decode("ascii")
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in [("coils.csv", coils(m)), ("vessel.csv", vessel(m)),
                       ("limiter.csv", limiter(eqdsk))]:
        (OUT / name).write_text("\n".join(HEADERS[name] + rows) + "\n", newline="\n")
        print(name, len(rows) - 1, "rows")


if __name__ == "__main__":
    main()
