#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace shafranov {

struct Coil {
    std::string name, circuit;
    double R, Z, dR, dZ, turns, resistance;
};
struct Filament {
    double R, Z, dR, dZ, resistance;
};
struct Point {
    double R, Z;
};
struct Machine {
    std::vector<Coil> coils;
    std::vector<Filament> vessel;
    std::vector<Point> limiter;
};

// Loads coils.csv, vessel.csv and limiter.csv from dir. Throws std::runtime_error naming the file
// (and the line, where there is one) on any error: missing file, wrong header, wrong field count,
// unparsable number, R <= 0, dR or dZ <= 0, resistance < 0, fewer than 3 limiter points,
// duplicate coil name.
Machine load_machine(const std::filesystem::path& dir);

} // namespace shafranov
