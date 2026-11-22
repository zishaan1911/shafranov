# Shafranov

A real-time, free-boundary tokamak plasma simulator for [TCV](https://www.epfl.ch/research/domains/swiss-plasma-center/research/tcv/),
the Tokamak à Configuration Variable at EPFL, written in C++20.

It solves the Grad–Shafranov equation on TCV's coil and vessel geometry and evolves the coil, vessel and plasma
currents in time, so the vertical instability comes out of the equations rather than being scripted.

**Status:** early. The project skeleton builds; the solver is next.

## Build

Needs CMake 3.24+, Ninja and a C++20 compiler (gcc 13+ tested).

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

## Validation

Every claim the simulator makes is backed by a test:

| Check | Status |
|---|---|
| Solov'ev analytic equilibrium, second-order convergence | planned |
| Green's functions against analytic limits | planned |
| Code-to-code against FreeGS on TCV geometry | planned |
| Vertical growth rates against published TCV values | planned |

## Licence

MIT
