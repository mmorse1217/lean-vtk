cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)

# Explicitly list the test source code and headers. The Catch header-only unit
# test framework is stored in with the test source.
set(LeanVTK_TEST_SRC
    test_lean_vtk.cpp
)
set(LeanVTK_TEST_HEADER
    catch.hpp
)

PREPEND(LeanVTK_TEST_SRC)

# Make an executable target that depends on the test source code we specified
# above.
add_executable(TestLeanVTK ${LeanVTK_TEST_SRC} ${LeanVTK_TEST_HEADER})

# Enable testing via CTest
enable_testing()
# Add our test as runnable via CTest
add_test(NAME TestLeanVTK  COMMAND TestLeanVTK)

# Link our unit tests against the library we compiled
target_link_libraries(TestLeanVTK LeanVTK)
