include(FetchContent)
# doctest's own CMakeLists predates CMake 3.5; CMake 4 needs this to configure it.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
FetchContent_Declare(doctest
  URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.12.tar.gz)
FetchContent_MakeAvailable(doctest)

add_executable(shafranov_tests
  test_main.cpp
  test_constants.cpp
  test_greens.cpp)
target_link_libraries(shafranov_tests PRIVATE shafranov_core shafranov_warnings doctest::doctest)
add_test(NAME shafranov_tests COMMAND shafranov_tests)
