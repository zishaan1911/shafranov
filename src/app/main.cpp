#include "app/view.h"
#include "core/circuits.h"
#include "core/equilibrium.h"
#include "core/grid.h"
#include "core/machine.h"
#include "sim/discharge.h"
#include "sim/worker.h"

#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <implot.h>

#include <cstdio>
#include <exception>
#include <filesystem>

int main(int argc, char** argv) {
    using namespace shafranov;

    const std::filesystem::path data_dir =
        argc > 1 ? std::filesystem::path(argv[1]) : std::filesystem::path(SHAFRANOV_DATA_DIR "/tcv");
    Machine machine;
    try {
        machine = load_machine(data_dir);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    const Grid grid(65, 129, 0.55, 1.20, -0.80, 0.80);
    const MachineGrid mg = make_machine_grid(machine, grid);
    const CircuitModel model = make_circuit_model(machine, mg);
    Discharge discharge(mg, model, demo_scenario());

    if (!glfwInit()) {
        return 1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* window = glfwCreateWindow(1400, 900, "Shafranov \xE2\x80\x94 TCV", nullptr, nullptr);
    if (window == nullptr) {
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGui::StyleColorsDark();
    ImPlot::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    {
        Worker worker(mg);
        ViewState state;
        worker.request(state.request);

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(viewport->WorkPos);
            ImGui::SetNextWindowSize(viewport->WorkSize);
            ImGui::Begin("Shafranov", nullptr,
                         ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                             ImGuiWindowFlags_NoSavedSettings |
                             ImGuiWindowFlags_NoBringToFrontOnFocus);
            draw(state, worker, discharge, mg, machine);
            ImGui::End();

            ImGui::Render();
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window);
        }
    } // the worker thread stops before the grid maps and ImGui go away

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwTerminate();
}
