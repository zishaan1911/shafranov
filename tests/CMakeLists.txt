include(FetchContent)
# doctest's own CMakeLists predates CMake 3.5; CMake 4 needs this to configure it.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
FetchContent_Declare(doctest
  URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.12.tar.gz)
FetchContent_MakeAvailable(doctest)

add_executable(shafranov_tests
  test_main.cpp
  test_constants.cpp
  test_coils.cpp
  test_contour.cpp
  test_greens.cpp
  test_spectral.cpp
  test_eqdsk.cpp
  test_machine.cpp
  test_tcv_data.cpp
  test_poisson.cpp
  test_profiles.cpp
  test_topology.cpp
  test_equilibrium.cpp
  test_inverse.cpp)
target_link_libraries(shafranov_tests PRIVATE shafranov_core shafranov_warnings doctest::doctest)
target_compile_definitions(shafranov_tests PRIVATE
  SHAFRANOV_DATA_DIR="${PROJECT_SOURCE_DIR}/data"
  SHAFRANOV_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
add_test(NAME shafranov_tests COMMAND shafranov_tests)
