cmake_minimum_required(VERSION 3.14)
project(LinearSolverProject VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing
enable_testing()

add_library(linear_solver_lib STATIC IMPORTED)
set_target_properties(linear_solver_lib PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/../LinearEquationSolver/cmake-build-debug/src/liblinear_solver_lib.a"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../LinearEquationSolver/src/"
)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

