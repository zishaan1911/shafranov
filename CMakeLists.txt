cmake_minimum_required(VERSION 3.24)
project(shafranov LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(shafranov_warnings INTERFACE)
target_compile_options(shafranov_warnings INTERFACE -Wall -Wextra -Wpedantic -Werror)

add_library(shafranov_core STATIC
  src/core/circuits.cpp
  src/core/coils.cpp
  src/core/contour.cpp
  src/core/dst.cpp
  src/core/eqdsk.cpp
  src/core/equilibrium.cpp
  src/core/greens.cpp
  src/core/grid.cpp
  src/core/inverse.cpp
  src/core/machine.cpp
  src/core/poisson.cpp
  src/core/profiles.cpp
  src/core/topology.cpp
  src/sim/worker.cpp)
target_include_directories(shafranov_core PUBLIC src)
include(FetchContent)
# Header-only; SOURCE_SUBDIR points nowhere so Eigen's own (large) CMake project is not added.
FetchContent_Declare(eigen
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
  SOURCE_SUBDIR none)
FetchContent_MakeAvailable(eigen)
target_include_directories(shafranov_core SYSTEM PUBLIC ${eigen_SOURCE_DIR})

find_package(Threads REQUIRED)
target_link_libraries(shafranov_core PUBLIC Threads::Threads PRIVATE shafranov_warnings)

option(SHAFRANOV_APP "Build the live window (needs OpenGL)" ON)
if(SHAFRANOV_APP)
  add_subdirectory(src/app)
endif()

enable_testing()
add_subdirectory(tests)
