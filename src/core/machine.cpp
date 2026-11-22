#include "core/machine.h"

#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace shafranov {

namespace {

struct Row {
    std::size_t line;
    std::vector<std::string> fields;
};

std::string_view trim(std::string_view s) {
    constexpr std::string_view blanks = " \t\r";
    const std::size_t first = s.find_first_not_of(blanks);
    if (first == std::string_view::npos) {
        return {};
    }
    return s.substr(first, s.find_last_not_of(blanks) - first + 1);
}

[[noreturn]] void fail(const std::filesystem::path& file, const std::string& what) {
    throw std::runtime_error(file.string() + ": " + what);
}

[[noreturn]] void fail(const std::filesystem::path& file, std::size_t line,
                       const std::string& what) {
    throw std::runtime_error(file.string() + " line " + std::to_string(line) + ": " + what);
}

std::vector<Row> read_rows(const std::filesystem::path& file, std::string_view header) {
    std::ifstream in(file);
    if (!in) {
        fail(file, "cannot open file");
    }
    std::vector<Row> rows;
    bool seen_header = false;
    std::string raw;
    std::size_t line = 0;
    while (std::getline(in, raw)) {
        ++line;
        const std::string_view text = trim(raw);
        if (text.empty() || text.front() == '#') {
            continue;
        }
        if (!seen_header) {
            if (text != header) {
                fail(file, line, "wrong header, expected '" + std::string(header) + "'");
            }
            seen_header = true;
            continue;
        }
        Row row{line, {}};
        std::size_t start = 0;
        while (true) {
            const std::size_t comma = text.find(',', start);
            if (comma == std::string_view::npos) {
                row.fields.emplace_back(trim(text.substr(start)));
                break;
            }
            row.fields.emplace_back(trim(text.substr(start, comma - start)));
            start = comma + 1;
        }
        rows.push_back(std::move(row));
    }
    if (!seen_header) {
        fail(file, "missing header");
    }
    return rows;
}

void expect_fields(const std::filesystem::path& file, const Row& row, std::size_t count) {
    if (row.fields.size() != count) {
        fail(file, row.line,
             "expected " + std::to_string(count) + " fields, got " +
                 std::to_string(row.fields.size()));
    }
}

double number(const std::filesystem::path& file, const Row& row, std::size_t index) {
    const std::string& s = row.fields[index];
    double value = 0.0;
    const char* end = s.data() + s.size();
    const auto [ptr, ec] = std::from_chars(s.data(), end, value);
    if (s.empty() || ec != std::errc{} || ptr != end || !std::isfinite(value)) {
        fail(file, row.line, "cannot parse number '" + s + "'");
    }
    return value;
}

// Positive-valued fields (R, dR, dZ) share one check so the messages stay uniform.
double positive(const std::filesystem::path& file, const Row& row, std::size_t index,
                const char* name) {
    const double value = number(file, row, index);
    if (value <= 0.0) {
        fail(file, row.line, std::string(name) + " must be > 0");
    }
    return value;
}

double resistance(const std::filesystem::path& file, const Row& row, std::size_t index) {
    const double value = number(file, row, index);
    if (value < 0.0) {
        fail(file, row.line, "resistance must be >= 0");
    }
    return value;
}

} // namespace

Machine load_machine(const std::filesystem::path& dir) {
    Machine machine;

    const auto coils_file = dir / "coils.csv";
    std::set<std::string> names;
    for (const Row& row : read_rows(coils_file, "name,circuit,R,Z,dR,dZ,turns,resistance")) {
        expect_fields(coils_file, row, 8);
        Coil coil;
        coil.name = row.fields[0];
        coil.circuit = row.fields[1];
        coil.R = positive(coils_file, row, 2, "R");
        coil.Z = number(coils_file, row, 3);
        coil.dR = positive(coils_file, row, 4, "dR");
        coil.dZ = positive(coils_file, row, 5, "dZ");
        coil.turns = number(coils_file, row, 6);
        coil.resistance = resistance(coils_file, row, 7);
        if (!names.insert(coil.name).second) {
            fail(coils_file, row.line, "duplicate coil name '" + coil.name + "'");
        }
        machine.coils.push_back(std::move(coil));
    }

    const auto vessel_file = dir / "vessel.csv";
    for (const Row& row : read_rows(vessel_file, "R,Z,dR,dZ,resistance")) {
        expect_fields(vessel_file, row, 5);
        Filament f;
        f.R = positive(vessel_file, row, 0, "R");
        f.Z = number(vessel_file, row, 1);
        f.dR = positive(vessel_file, row, 2, "dR");
        f.dZ = positive(vessel_file, row, 3, "dZ");
        f.resistance = resistance(vessel_file, row, 4);
        machine.vessel.push_back(f);
    }

    const auto limiter_file = dir / "limiter.csv";
    for (const Row& row : read_rows(limiter_file, "R,Z")) {
        expect_fields(limiter_file, row, 2);
        Point p;
        p.R = positive(limiter_file, row, 0, "R");
        p.Z = number(limiter_file, row, 1);
        machine.limiter.push_back(p);
    }
    if (machine.limiter.size() < 3) {
        fail(limiter_file,
             "need at least 3 points, got " + std::to_string(machine.limiter.size()));
    }

    return machine;
}

} // namespace shafranov
