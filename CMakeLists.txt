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
  src/core/greens.cpp
  src/core/grid.cpp)
target_include_directories(shafranov_core PUBLIC src)
target_link_libraries(shafranov_core PRIVATE shafranov_warnings)

enable_testing()
add_subdirectory(tests)
