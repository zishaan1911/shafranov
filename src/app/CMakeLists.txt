include(FetchContent)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_Declare(glfw URL https://github.com/glfw/glfw/archive/refs/tags/3.4.tar.gz)
FetchContent_Declare(imgui URL https://github.com/ocornut/imgui/archive/refs/tags/v1.92.9.tar.gz)
FetchContent_Declare(implot URL https://github.com/epezent/implot/archive/refs/tags/v1.1.tar.gz)
FetchContent_MakeAvailable(glfw imgui implot)
find_package(OpenGL REQUIRED)

# Dear ImGui and ImPlot ship no CMake targets; build them as one library.
add_library(imgui_implot STATIC
  ${imgui_SOURCE_DIR}/imgui.cpp
  ${imgui_SOURCE_DIR}/imgui_draw.cpp
  ${imgui_SOURCE_DIR}/imgui_tables.cpp
  ${imgui_SOURCE_DIR}/imgui_widgets.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
  ${implot_SOURCE_DIR}/implot.cpp
  ${implot_SOURCE_DIR}/implot_items.cpp)
target_include_directories(imgui_implot SYSTEM PUBLIC
  ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends ${implot_SOURCE_DIR})
target_link_libraries(imgui_implot PUBLIC glfw OpenGL::GL)

add_executable(shafranov_app main.cpp)
target_link_libraries(shafranov_app PRIVATE shafranov_core imgui_implot shafranov_warnings)
if(MINGW)
  # Runs on double-click without the MSYS2 runtime DLLs on PATH.
  target_link_options(shafranov_app PRIVATE -static)
endif()
