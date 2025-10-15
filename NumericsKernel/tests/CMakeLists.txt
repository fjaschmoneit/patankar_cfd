cmake_minimum_required(VERSION 3.14)
project(NumeriksKernel_tests VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing
enable_testing()

# Download and configure GoogleTest
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

# Create test executable
add_executable(linear_solver_tests
        test_linearEquationSolver.cpp
        test_performance.cpp
)

# Link with GoogleTest and our static library
target_link_libraries(linear_solver_tests
        PRIVATE
        gtest_main
        linear_solver_lib
)

include_directories(/usr/local/include)

# Include GoogleTest
include(GoogleTest)
gtest_discover_tests(linear_solver_tests)
