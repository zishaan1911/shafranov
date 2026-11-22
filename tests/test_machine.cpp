#include "core/machine.h"

#include <doctest/doctest.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

namespace fs = std::filesystem;

const std::string valid_coils =
    "# poloidal field coils\r\n"
    "name,circuit,R,Z,dR,dZ,turns,resistance\r\n"
    "\r\n"
    "E1,E,1.5,0.5,0.1,0.2,10,0.01\r\n"
    "G1, G ,2.0,-1.0,0.2,0.2,5.5,0.02\r\n"
    "# anti-series\r\n"
    "G2,G,2.0,1.0,0.2,0.2,-5.5,0.02\r\n";
const std::string valid_vessel =
    "R,Z,dR,dZ,resistance\r\n"
    "1.8,0.0,0.05,0.5,1e-4\r\n"
    "# lower\r\n"
    "1.8,-0.5,0.05,0.5,2e-4\r\n";
const std::string valid_limiter =
    "R,Z\r\n"
    "1.0,-1.0\r\n"
    "2.0,-1.0\r\n"
    "2.0,1.0\r\n"
    "1.0,1.0\r\n";

// Owns a fresh directory for one test case and removes it at the end.
struct TempDir {
    fs::path path;

    explicit TempDir(const std::string& name)
        : path(fs::temp_directory_path() / ("shafranov_test_machine_" + name)) {
        fs::remove_all(path);
        fs::create_directories(path);
    }
    ~TempDir() {
        std::error_code ec;
        fs::remove_all(path, ec);
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;

    void write(const std::string& file, const std::string& text) const {
        std::ofstream out(path / file, std::ios::binary);
        out << text;
    }
    void write_valid() const {
        write("coils.csv", valid_coils);
        write("vessel.csv", valid_vessel);
        write("limiter.csv", valid_limiter);
    }
};

std::string load_error(const fs::path& dir) {
    try {
        shafranov::load_machine(dir);
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "<no error>";
}

bool contains(const std::string& text, const std::string& part) {
    return text.find(part) != std::string::npos;
}

} // namespace

TEST_CASE("load_machine reads a valid machine") {
    TempDir dir("valid");
    dir.write_valid();
    const auto m = shafranov::load_machine(dir.path);

    REQUIRE(m.coils.size() == 3);
    CHECK(m.coils[0].name == "E1");
    CHECK(m.coils[0].circuit == "E");
    CHECK(m.coils[0].R == doctest::Approx(1.5));
    CHECK(m.coils[0].Z == doctest::Approx(0.5));
    CHECK(m.coils[0].dR == doctest::Approx(0.1));
    CHECK(m.coils[0].dZ == doctest::Approx(0.2));
    CHECK(m.coils[0].turns == doctest::Approx(10.0));
    CHECK(m.coils[0].resistance == doctest::Approx(0.01));
    CHECK(m.coils[1].name == "G1");
    CHECK(m.coils[1].circuit == "G");
    CHECK(m.coils[1].Z == doctest::Approx(-1.0));
    CHECK(m.coils[1].turns == doctest::Approx(5.5));
    CHECK(m.coils[2].name == "G2");
    CHECK(m.coils[2].circuit == "G");
    CHECK(m.coils[2].Z == doctest::Approx(1.0));
    CHECK(m.coils[2].turns == doctest::Approx(-5.5));
    CHECK(m.coils[2].resistance == doctest::Approx(0.02));

    REQUIRE(m.vessel.size() == 2);
    CHECK(m.vessel[0].R == doctest::Approx(1.8));
    CHECK(m.vessel[0].Z == doctest::Approx(0.0));
    CHECK(m.vessel[0].dR == doctest::Approx(0.05));
    CHECK(m.vessel[0].dZ == doctest::Approx(0.5));
    CHECK(m.vessel[0].resistance == doctest::Approx(1e-4));
    CHECK(m.vessel[1].Z == doctest::Approx(-0.5));
    CHECK(m.vessel[1].resistance == doctest::Approx(2e-4));

    REQUIRE(m.limiter.size() == 4);
    CHECK(m.limiter[0].R == doctest::Approx(1.0));
    CHECK(m.limiter[0].Z == doctest::Approx(-1.0));
    CHECK(m.limiter[1].R == doctest::Approx(2.0));
    CHECK(m.limiter[2].Z == doctest::Approx(1.0));
    CHECK(m.limiter[3].R == doctest::Approx(1.0));
    CHECK(m.limiter[3].Z == doctest::Approx(1.0));
}

TEST_CASE("load_machine rejects a wrong header") {
    TempDir dir("header");
    dir.write_valid();
    dir.write("coils.csv", "name,circuit,R,Z\n");
    CHECK(contains(load_error(dir.path), "coils.csv"));
}

TEST_CASE("load_machine rejects an unparsable number") {
    TempDir dir("number");
    dir.write_valid();
    dir.write("vessel.csv",
              "# comment\n"
              "R,Z,dR,dZ,resistance\n"
              "1.8,0.0,0.05,0.5,1e-4\n"
              "1.8,abc,0.05,0.5,2e-4\n");
    const std::string msg = load_error(dir.path);
    CHECK(contains(msg, "vessel.csv"));
    CHECK(contains(msg, "line 4"));
}

TEST_CASE("load_machine rejects trailing junk after a number") {
    TempDir dir("junk");
    dir.write_valid();
    dir.write("limiter.csv", "R,Z\n1.0,1.0x\n2.0,1.0\n2.0,2.0\n");
    const std::string msg = load_error(dir.path);
    CHECK(contains(msg, "limiter.csv"));
    CHECK(contains(msg, "line 2"));
}

TEST_CASE("load_machine rejects R <= 0") {
    TempDir dir("radius");
    dir.write_valid();
    dir.write("coils.csv",
              "name,circuit,R,Z,dR,dZ,turns,resistance\n"
              "E1,E,0.0,0.5,0.1,0.2,10,0.01\n");
    CHECK(contains(load_error(dir.path), "coils.csv"));
}

TEST_CASE("load_machine rejects a duplicate coil name") {
    TempDir dir("duplicate");
    dir.write_valid();
    dir.write("coils.csv",
              "name,circuit,R,Z,dR,dZ,turns,resistance\n"
              "E1,E,1.5,0.5,0.1,0.2,10,0.01\n"
              "E1,E,1.6,0.5,0.1,0.2,10,0.01\n");
    CHECK(contains(load_error(dir.path), "coils.csv"));
}

TEST_CASE("load_machine rejects a missing vessel.csv") {
    TempDir dir("missing");
    dir.write_valid();
    fs::remove(dir.path / "vessel.csv");
    CHECK(contains(load_error(dir.path), "vessel.csv"));
}

TEST_CASE("load_machine rejects fewer than 3 limiter points") {
    TempDir dir("limiter");
    dir.write_valid();
    dir.write("limiter.csv", "R,Z\n1.0,-1.0\n2.0,-1.0\n");
    CHECK(contains(load_error(dir.path), "limiter.csv"));
}
