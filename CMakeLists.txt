cmake_minimum_required(VERSION 3.14)
project(Patankar_CFD VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing, I only need this, when I am executing test in cmake, But I can alo simply use CLion's test environment
# enable_testing()

add_subdirectory(FVM)
add_subdirectory(LinearEquationSolver)
