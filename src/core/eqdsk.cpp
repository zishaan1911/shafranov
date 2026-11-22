#include "core/eqdsk.h"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace shafranov {

namespace {

constexpr std::size_t field_width = 16;
constexpr std::size_t fields_per_line = 5;
constexpr std::size_t description_width = 48;

bool is_blank(const std::string& s) {
    return s.find_first_not_of(" \t\r\n") == std::string::npos;
}

// ---- writing ----

std::string format_real(double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%16.9E", v);
    return buf;
}

void write_values(std::ostream& out, const std::vector<double>& values) {
    for (std::size_t k = 0; k < values.size(); ++k) {
        out << format_real(values[k]);
        if ((k + 1) % fields_per_line == 0 || k + 1 == values.size())
            out << '\n';
    }
}

std::vector<double> interleave(const std::vector<double>& r, const std::vector<double>& z) {
    std::vector<double> pairs;
    pairs.reserve(2 * r.size());
    for (std::size_t k = 0; k < r.size(); ++k) {
        pairs.push_back(r[k]);
        pairs.push_back(z[k]);
    }
    return pairs;
}

void require(bool ok, const char* message) {
    if (!ok)
        throw std::invalid_argument(message);
}

// ---- reading ----

struct LineReader {
    std::istream& in;

    std::string next_line(const char* what) {
        std::string line;
        if (!std::getline(in, line))
            throw std::runtime_error(std::string("G-EQDSK: unexpected end of file while reading ") +
                                     what);
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return line;
    }

    // Reads `count` reals from fresh lines, parsing fixed 16-character fields.
    std::vector<double> values(std::size_t count, const char* what) {
        std::vector<double> out;
        while (out.size() < count) {
            const std::string line = next_line(what);
            for (std::size_t pos = 0; pos < line.size() && out.size() < count;
                 pos += field_width) {
                const std::string field = line.substr(pos, field_width);
                if (is_blank(field))
                    continue;
                char* end = nullptr;
                const double v = std::strtod(field.c_str(), &end);
                if (end == field.c_str() || !is_blank(end))
                    throw std::runtime_error(std::string("G-EQDSK: bad number '") + field +
                                             "' while reading " + what);
                out.push_back(v);
            }
        }
        return out;
    }
};

int parse_int(const std::string& text, const char* what) {
    char* end = nullptr;
    const long v = std::strtol(text.c_str(), &end, 10);
    if (end == text.c_str() || !is_blank(end))
        throw std::runtime_error(std::string("G-EQDSK: bad integer '") + text + "' in " + what);
    return static_cast<int>(v);
}

std::string rtrim(std::string s) {
    s.erase(s.find_last_not_of(' ') + 1);
    return s;
}

} // namespace

void write_geqdsk(std::ostream& out, const GEqdsk& eq) {
    require(eq.nw >= 0 && eq.nh >= 0, "write_geqdsk: nw and nh must be non-negative");
    const std::size_t nw = static_cast<std::size_t>(eq.nw);
    const std::size_t nh = static_cast<std::size_t>(eq.nh);
    require(eq.fpol.size() == nw, "write_geqdsk: fpol.size() != nw");
    require(eq.pres.size() == nw, "write_geqdsk: pres.size() != nw");
    require(eq.ffprim.size() == nw, "write_geqdsk: ffprim.size() != nw");
    require(eq.pprime.size() == nw, "write_geqdsk: pprime.size() != nw");
    require(eq.qpsi.size() == nw, "write_geqdsk: qpsi.size() != nw");
    require(eq.psirz.size() == nw * nh, "write_geqdsk: psirz.size() != nw * nh");
    require(eq.rbbbs.size() == eq.zbbbs.size(), "write_geqdsk: rbbbs.size() != zbbbs.size()");
    require(eq.rlim.size() == eq.zlim.size(), "write_geqdsk: rlim.size() != zlim.size()");

    std::string desc = eq.description.substr(0, description_width);
    desc.resize(description_width, ' ');
    char head[32];
    std::snprintf(head, sizeof head, "%4d%4d%4d", 0, eq.nw, eq.nh);
    out << desc << head << '\n';

    write_values(out, {eq.rdim, eq.zdim, eq.rcentr, eq.rleft, eq.zmid});
    write_values(out, {eq.rmaxis, eq.zmaxis, eq.simag, eq.sibry, eq.bcentr});
    write_values(out, {eq.current, eq.simag, 0.0, eq.rmaxis, 0.0});
    write_values(out, {eq.zmaxis, 0.0, eq.sibry, 0.0, 0.0});

    write_values(out, eq.fpol);
    write_values(out, eq.pres);
    write_values(out, eq.ffprim);
    write_values(out, eq.pprime);

    // File order is R fastest; memory order is i * nh + j.
    std::vector<double> file_order;
    file_order.reserve(eq.psirz.size());
    for (std::size_t j = 0; j < nh; ++j)
        for (std::size_t i = 0; i < nw; ++i)
            file_order.push_back(eq.psirz[i * nh + j]);
    write_values(out, file_order);

    write_values(out, eq.qpsi);

    char counts[32];
    std::snprintf(counts, sizeof counts, "%5d%5d", static_cast<int>(eq.rbbbs.size()),
                  static_cast<int>(eq.rlim.size()));
    out << counts << '\n';
    write_values(out, interleave(eq.rbbbs, eq.zbbbs));
    write_values(out, interleave(eq.rlim, eq.zlim));
}

GEqdsk read_geqdsk(std::istream& in) {
    LineReader reader{in};
    GEqdsk eq;

    const std::string header = reader.next_line("header");
    eq.description = rtrim(header.substr(0, description_width));
    std::istringstream tokens(header);
    std::vector<std::string> words;
    for (std::string w; tokens >> w;)
        words.push_back(w);
    if (words.size() < 2)
        throw std::runtime_error("G-EQDSK: header line has no nw and nh");
    eq.nw = parse_int(words[words.size() - 2], "header (nw)");
    eq.nh = parse_int(words[words.size() - 1], "header (nh)");
    if (eq.nw <= 0 || eq.nh <= 0)
        throw std::runtime_error("G-EQDSK: nw and nh must be positive");
    const std::size_t nw = static_cast<std::size_t>(eq.nw);
    const std::size_t nh = static_cast<std::size_t>(eq.nh);

    const auto l1 = reader.values(fields_per_line, "scalar line 1");
    eq.rdim = l1[0];
    eq.zdim = l1[1];
    eq.rcentr = l1[2];
    eq.rleft = l1[3];
    eq.zmid = l1[4];
    const auto l2 = reader.values(fields_per_line, "scalar line 2");
    eq.rmaxis = l2[0];
    eq.zmaxis = l2[1];
    eq.simag = l2[2];
    eq.sibry = l2[3];
    eq.bcentr = l2[4];
    const auto l3 = reader.values(fields_per_line, "scalar line 3");
    eq.current = l3[0];
    reader.values(fields_per_line, "scalar line 4");

    eq.fpol = reader.values(nw, "fpol");
    eq.pres = reader.values(nw, "pres");
    eq.ffprim = reader.values(nw, "ffprim");
    eq.pprime = reader.values(nw, "pprime");

    const auto file_order = reader.values(nw * nh, "psirz");
    eq.psirz.resize(nw * nh);
    for (std::size_t j = 0; j < nh; ++j)
        for (std::size_t i = 0; i < nw; ++i)
            eq.psirz[i * nh + j] = file_order[j * nw + i];

    eq.qpsi = reader.values(nw, "qpsi");

    const std::string counts = reader.next_line("boundary counts");
    if (counts.size() <= 5)
        throw std::runtime_error("G-EQDSK: boundary count line must hold nbbbs and limitr");
    const int nbbbs = parse_int(counts.substr(0, 5), "nbbbs");
    const int limitr = parse_int(counts.substr(5, 5), "limitr");
    if (nbbbs < 0 || limitr < 0)
        throw std::runtime_error("G-EQDSK: nbbbs and limitr must be non-negative");

    const auto split = [](const std::vector<double>& pairs, std::vector<double>& r,
                          std::vector<double>& z) {
        for (std::size_t k = 0; k + 1 < pairs.size(); k += 2) {
            r.push_back(pairs[k]);
            z.push_back(pairs[k + 1]);
        }
    };
    split(reader.values(2 * static_cast<std::size_t>(nbbbs), "boundary"), eq.rbbbs, eq.zbbbs);
    split(reader.values(2 * static_cast<std::size_t>(limitr), "limiter"), eq.rlim, eq.zlim);
    return eq;
}

} // namespace shafranov
