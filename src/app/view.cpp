#include "app/view.h"

#include "core/contour.h"

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>

namespace shafranov {

namespace {

constexpr double supply_limit_kA = 7.7;

bool is_ef(const std::string& name) {
    return name.size() >= 2 && (name[0] == 'E' || name[0] == 'F') && name[1] >= '0' &&
           name[1] <= '9';
}

// ImPlot item flags and the per-plot flags are separate unscoped enums; combine as ints.
int flags(int a, int b = 0) { return a | b; }

// Slider on a double shown as value * scale.
bool slider(const char* label, double& v, double scale, float lo, float hi, const char* fmt) {
    float f = static_cast<float>(v * scale);
    const bool changed = ImGui::SliderFloat(label, &f, lo, hi, fmt);
    if (changed) {
        v = static_cast<double>(f) / scale;
    }
    return changed;
}

void add_segments(std::vector<double>& xs, std::vector<double>& ys, const std::vector<Segment>& s) {
    for (const Segment& seg : s) {
        xs.push_back(seg.R1);
        ys.push_back(seg.Z1);
        xs.push_back(seg.R2);
        ys.push_back(seg.Z2);
    }
}

void rebuild_caches(ViewState& st, const Snapshot& snap, const MachineGrid& mg) {
    st.heat.clear();
    st.flux_x.clear();
    st.flux_y.clear();
    st.sep_x.clear();
    st.sep_y.clear();
    const Grid& g = mg.grid;
    const Equilibrium& eq = snap.eq;
    if (!snap.error.empty() || eq.psi.size() != static_cast<std::size_t>(g.size())) {
        return;
    }
    const double psi_axis = eq.topology.axis.psi;
    const double span = eq.topology.psi_bndry - psi_axis;
    if (std::abs(span) < 1e-30) {
        return;
    }

    std::vector<double> pn(eq.psi.size());
    for (std::size_t k = 0; k < pn.size(); ++k) {
        pn[k] = (eq.psi[k] - psi_axis) / span;
    }

    st.heat.resize(pn.size());
    for (int j = 0; j < g.nZ; ++j) {
        for (int i = 0; i < g.nR; ++i) {
            const double v = std::clamp(pn[static_cast<std::size_t>(g.idx(i, j))], 0.0, 2.0);
            st.heat[static_cast<std::size_t>((g.nZ - 1 - j) * g.nR + i)] = static_cast<float>(v);
        }
    }

    for (int k = 1; k <= 9; ++k) {
        add_segments(st.flux_x, st.flux_y, contour(g, pn, 0.1 * k));
    }
    for (const double level : {1.05, 1.1, 1.2, 1.4}) {
        add_segments(st.flux_x, st.flux_y, contour(g, pn, level));
    }
    add_segments(st.sep_x, st.sep_y, contour(g, pn, 1.0));
}

ImU32 coil_color(double current_kA) {
    const double t = std::clamp(current_kA / supply_limit_kA, -1.0, 1.0);
    const double grey = 110.0;
    const double m = std::abs(t);
    const double r = grey + m * ((t > 0 ? 220.0 : 40.0) - grey);
    const double gr = grey + m * (50.0 - grey);
    const double b = grey + m * ((t > 0 ? 40.0 : 220.0) - grey);
    return IM_COL32(static_cast<int>(r), static_cast<int>(gr), static_cast<int>(b), 255);
}

void draw_controls(ViewState& st, Worker& worker, const MachineGrid& mg,
                   const std::shared_ptr<const Snapshot>& snap) {
    Request& r = st.request;
    bool changed = false;

    ImGui::SeparatorText("Plasma");
    changed |= slider("Ip (kA)", r.Ip, 1e-3, 50.0f, 400.0f, "%.0f");
    changed |= slider("beta0", r.beta0, 1.0, 0.0f, 1.0f, "%.2f");
    changed |= slider("R0 (m)", r.R0, 1.0, 0.80f, 0.95f, "%.3f");
    changed |= slider("a (m)", r.a, 1.0, 0.12f, 0.25f, "%.3f");
    changed |= slider("kappa", r.kappa, 1.0, 1.0f, 2.0f, "%.2f");
    changed |= slider("delta", r.delta, 1.0, -0.5f, 0.6f, "%.2f");
    changed |= ImGui::Checkbox("Diverted (lower X-point)", &r.diverted);
    if (r.diverted) {
        changed |= slider("X-point R (m)", r.xR, 1.0, 0.70f, 1.00f, "%.3f");
        changed |= slider("X-point Z (m)", r.xZ, 1.0, -0.60f, -0.25f, "%.3f");
    }

    ImGui::SeparatorText("Coils");
    if (st.slider.size() != mg.circuits.size()) {
        st.slider.assign(mg.circuits.size(), 0.0f);
    }
    int mode = st.fixed_currents ? 1 : 0;
    bool clicked = ImGui::RadioButton("Shape control", &mode, 0);
    ImGui::SameLine();
    clicked |= ImGui::RadioButton("Fixed currents", &mode, 1);
    if (clicked && (mode == 1) != st.fixed_currents) {
        st.fixed_currents = (mode == 1);
        if (st.fixed_currents && snap && snap->currents.circuits.size() == mg.circuits.size()) {
            for (std::size_t k = 0; k < mg.circuits.size(); ++k) {
                st.slider[k] = static_cast<float>(snap->currents.circuits[k] * 1e-3);
            }
        }
        changed = true;
    }
    if (st.fixed_currents) {
        for (std::size_t k = 0; k < mg.circuits.size(); ++k) {
            if (is_ef(mg.circuits[k])) {
                const float lim = static_cast<float>(supply_limit_kA);
                changed |= ImGui::SliderFloat(mg.circuits[k].c_str(), &st.slider[k], -lim, lim,
                                              "%.2f kA");
            }
        }
        ImGui::TextWrapped("Fixed currents give no shape or vertical control: elongated plasmas "
                           "will drift and may be lost.");
    }

    if (changed) {
        Request req = r;
        if (st.fixed_currents) {
            // Circuits without a slider (OH, G, T) keep the latest solution's value.
            std::vector<double> c(mg.circuits.size(), 0.0);
            if (snap && snap->currents.circuits.size() == c.size()) {
                c = snap->currents.circuits;
            }
            for (std::size_t k = 0; k < c.size(); ++k) {
                if (is_ef(mg.circuits[k])) {
                    c[k] = static_cast<double>(st.slider[k]) * 1e3;
                }
            }
            req.fixed_currents = c;
        } else {
            req.fixed_currents.reset();
        }
        worker.request(req);
    }

    ImGui::SeparatorText("Solution");
    if (!snap) {
        ImGui::TextUnformatted("no solution yet");
    } else if (!snap->error.empty()) {
        ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f), "%s", snap->error.c_str());
    } else if (snap->eq.converged) {
        ImGui::TextColored(ImVec4(0.3f, 0.9f, 0.3f, 1.0f), "converged");
    } else {
        ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.2f, 1.0f), "did not converge");
    }
    if (worker.busy()) {
        ImGui::TextUnformatted("solving...");
    }
    if (snap && snap->error.empty()) {
        const Topology& t = snap->eq.topology;
        ImGui::Text("iterations: %d", snap->eq.iterations);
        ImGui::Text("solve time: %.1f ms", snap->solve_ms);
        ImGui::Text("axis: R = %.3f m, Z = %.3f m", t.axis.R, t.axis.Z);
        ImGui::Text("psi_axis: %.4f Wb/rad", t.axis.psi);
        ImGui::Text("psi_bndry: %.4f Wb/rad", t.psi_bndry);
        ImGui::TextUnformatted(t.limited ? "limited" : "diverted");
        if (snap->metrics) {
            const ShapeMetrics& m = *snap->metrics;
            ImGui::Text("kappa: %.3f", m.kappa);
            ImGui::Text("delta: %.3f upper, %.3f lower", m.delta_upper, m.delta_lower);
            ImGui::Text("a: %.3f m, Rgeo: %.3f m", m.a, m.Rgeo);
        }
    }

    ImGui::SeparatorText("Coil currents");
    std::vector<std::size_t> ef;
    for (std::size_t k = 0; k < mg.circuits.size(); ++k) {
        if (is_ef(mg.circuits[k])) {
            ef.push_back(k);
        }
    }
    if (ImPlot::BeginPlot("##coil_currents", ImVec2(-1, 220), ImPlotFlags_NoLegend)) {
        std::vector<double> ticks, ok_x, ok_y, over_x, over_y;
        std::vector<const char*> labels;
        const bool have = snap && snap->currents.circuits.size() == mg.circuits.size();
        for (std::size_t n = 0; n < ef.size(); ++n) {
            const double x = static_cast<double>(n);
            ticks.push_back(x);
            labels.push_back(mg.circuits[ef[n]].c_str());
            if (!have) {
                continue;
            }
            const double kA = snap->currents.circuits[ef[n]] * 1e-3;
            auto& xs = std::abs(kA) > supply_limit_kA ? over_x : ok_x;
            auto& ys = std::abs(kA) > supply_limit_kA ? over_y : ok_y;
            xs.push_back(x);
            ys.push_back(kA);
        }
        ImPlot::SetupAxes(nullptr, "I (kA)", 0, 0);
        ImPlot::SetupAxisLimits(ImAxis_X1, -0.7, static_cast<double>(ef.size()) - 0.3,
                                ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, -9.0, 9.0, ImPlotCond_Always);
        ImPlot::SetupAxisTicks(ImAxis_X1, ticks.data(), static_cast<int>(ticks.size()),
                               labels.data());
        if (!ok_x.empty()) {
            ImPlotSpec spec;
            spec.FillColor = ImVec4(0.30f, 0.55f, 0.90f, 1.0f);
            ImPlot::PlotBars("##ok", ok_x.data(), ok_y.data(), static_cast<int>(ok_x.size()), 0.67,
                             spec);
        }
        if (!over_x.empty()) {
            ImPlotSpec spec;
            spec.FillColor = ImVec4(0.90f, 0.20f, 0.20f, 1.0f);
            ImPlot::PlotBars("##over", over_x.data(), over_y.data(),
                             static_cast<int>(over_x.size()), 0.67, spec);
        }
        const double limits[2] = {supply_limit_kA, -supply_limit_kA};
        ImPlotSpec lim;
        lim.LineColor = ImVec4(0.95f, 0.75f, 0.2f, 1.0f);
        lim.Flags = flags(ImPlotItemFlags_NoLegend, ImPlotInfLinesFlags_Horizontal);
        ImPlot::PlotInfLines("##limit", limits, 2, lim);
        ImPlot::EndPlot();
    }
}

void draw_plot(ViewState& st, const MachineGrid& mg, const Machine& machine,
               const std::shared_ptr<const Snapshot>& snap) {
    const Grid& g = mg.grid;
    if (!ImPlot::BeginPlot("##section", ImVec2(-1, -1),
                           flags(ImPlotFlags_Equal, ImPlotFlags_NoLegend))) {
        return;
    }
    ImPlot::SetupAxes("R (m)", "Z (m)", 0, 0);
    ImPlot::SetupAxesLimits(0.35, 1.45, -0.95, 0.95, ImPlotCond_Always);

    const bool have_field = snap && !st.heat.empty();
    if (have_field) {
        ImPlot::PushColormap(ImPlotColormap_Viridis);
        ImPlot::PlotHeatmap("##psi_n", st.heat.data(), g.nZ, g.nR, 0.0, 2.0, nullptr,
                            ImPlotPoint(g.Rmin, g.Zmin), ImPlotPoint(g.Rmax, g.Zmax));
        ImPlot::PopColormap();

        const int no_legend = ImPlotItemFlags_NoLegend;
        if (!st.flux_x.empty()) {
            ImPlotSpec spec;
            spec.LineColor = ImVec4(1.0f, 1.0f, 1.0f, 0.45f);
            spec.LineWeight = 1.0f;
            spec.Flags = flags(no_legend, ImPlotLineFlags_Segments);
            ImPlot::PlotLine("##flux", st.flux_x.data(), st.flux_y.data(),
                             static_cast<int>(st.flux_x.size()), spec);
        }
        if (!st.sep_x.empty()) {
            ImPlotSpec spec;
            spec.LineColor = ImVec4(1.0f, 1.0f, 0.6f, 1.0f);
            spec.LineWeight = 3.0f;
            spec.Flags = flags(no_legend, ImPlotLineFlags_Segments);
            ImPlot::PlotLine("##sep", st.sep_x.data(), st.sep_y.data(),
                             static_cast<int>(st.sep_x.size()), spec);
        }
    }

    if (!machine.limiter.empty()) {
        std::vector<double> lx, ly;
        for (const Point& p : machine.limiter) {
            lx.push_back(p.R);
            ly.push_back(p.Z);
        }
        lx.push_back(machine.limiter.front().R);
        ly.push_back(machine.limiter.front().Z);
        ImPlotSpec spec;
        spec.LineColor = ImVec4(0.6f, 0.6f, 0.6f, 1.0f);
        spec.LineWeight = 2.5f;
        spec.Flags = ImPlotItemFlags_NoLegend;
        ImPlot::PlotLine("##limiter", lx.data(), ly.data(), static_cast<int>(lx.size()), spec);
    }

    if (have_field) {
        const Topology& t = snap->eq.topology;
        const int no_legend = ImPlotItemFlags_NoLegend;
        const double ax = t.axis.R, az = t.axis.Z;
        ImPlotSpec axis_spec;
        axis_spec.Marker = ImPlotMarker_Plus;
        axis_spec.MarkerSize = 7.0f;
        axis_spec.MarkerLineColor = ImVec4(1.0f, 1.0f, 1.0f, 1.0f);
        axis_spec.LineWeight = 2.0f;
        axis_spec.Flags = no_legend;
        ImPlot::PlotScatter("##axis", &ax, &az, 1, axis_spec);

        if (!t.xpoints.empty()) {
            std::vector<double> xr, xz;
            for (const CriticalPoint& p : t.xpoints) {
                xr.push_back(p.R);
                xz.push_back(p.Z);
            }
            ImPlotSpec spec;
            spec.Marker = ImPlotMarker_Cross;
            spec.MarkerSize = 7.0f;
            spec.MarkerLineColor = ImVec4(1.0f, 0.4f, 0.4f, 1.0f);
            spec.LineWeight = 2.0f;
            spec.Flags = no_legend;
            ImPlot::PlotScatter("##xpoints", xr.data(), xz.data(), static_cast<int>(xr.size()),
                                spec);
        }

        const ShapeTargets targets = targets_for(snap->request);
        std::vector<double> tr, tz;
        for (const auto& iso : targets.isoflux) {
            tr.push_back(iso.R1);
            tz.push_back(iso.Z1);
            tr.push_back(iso.R2);
            tz.push_back(iso.Z2);
        }
        for (const Point& p : targets.xpoints) {
            tr.push_back(p.R);
            tz.push_back(p.Z);
        }
        if (!tr.empty()) {
            ImPlotSpec spec;
            spec.Marker = ImPlotMarker_Circle;
            spec.MarkerSize = 4.0f;
            spec.MarkerLineColor = ImVec4(1.0f, 0.6f, 0.1f, 1.0f);
            spec.FillAlpha = 0.0f; // hollow
            spec.Flags = no_legend;
            ImPlot::PlotScatter("##targets", tr.data(), tz.data(), static_cast<int>(tr.size()),
                                spec);
        }
    }

    // Vessel, coils and the error message go straight on the plot's draw list.
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlot::PushPlotClipRect();
    for (const Filament& f : machine.vessel) {
        const ImVec2 lo = ImPlot::PlotToPixels(f.R - 0.5 * f.dR, f.Z + 0.5 * f.dZ);
        const ImVec2 hi = ImPlot::PlotToPixels(f.R + 0.5 * f.dR, f.Z - 0.5 * f.dZ);
        dl->AddRectFilled(lo, hi, IM_COL32(150, 150, 150, 255));
    }
    const bool have_currents = snap && snap->currents.circuits.size() == mg.circuits.size();
    for (const Coil& c : machine.coils) {
        double kA = 0.0;
        if (have_currents) {
            const auto it = std::find(mg.circuits.begin(), mg.circuits.end(), c.circuit);
            if (it != mg.circuits.end()) {
                kA = snap->currents.circuits[static_cast<std::size_t>(it - mg.circuits.begin())] *
                     1e-3;
            }
        }
        const ImVec2 lo = ImPlot::PlotToPixels(c.R - 0.5 * c.dR, c.Z + 0.5 * c.dZ);
        const ImVec2 hi = ImPlot::PlotToPixels(c.R + 0.5 * c.dR, c.Z - 0.5 * c.dZ);
        dl->AddRectFilled(lo, hi, coil_color(kA));
        dl->AddRect(lo, hi, IM_COL32(255, 255, 255, 255));
        if (is_ef(c.circuit)) {
            dl->AddText(ImVec2(hi.x + 3.0f, lo.y), IM_COL32(255, 255, 255, 255), c.name.c_str());
        }
    }
    if (snap && !snap->error.empty()) {
        const ImVec2 pos = ImPlot::GetPlotPos();
        dl->AddText(ImVec2(pos.x + 10.0f, pos.y + 8.0f), IM_COL32(255, 70, 70, 255),
                    snap->error.c_str());
    }
    ImPlot::PopPlotClipRect();

    ImPlot::EndPlot();
}

} // namespace

void draw(ViewState& state, Worker& worker, const MachineGrid& mg, const Machine& machine) {
    const std::shared_ptr<const Snapshot> snap = worker.latest();
    if (snap && snap->id != state.seen_id) {
        rebuild_caches(state, *snap, mg);
        state.seen_id = snap->id;
    }

    if (ImGui::BeginChild("controls", ImVec2(360.0f, 0.0f), ImGuiChildFlags_Borders)) {
        draw_controls(state, worker, mg, snap);
    }
    ImGui::EndChild();
    ImGui::SameLine();
    draw_plot(state, mg, machine, snap);
}

} // namespace shafranov
