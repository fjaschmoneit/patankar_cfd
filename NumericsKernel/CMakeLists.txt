cmake_minimum_required(VERSION 3.14)
project(numerikKernel VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)
#
#add_custom_target(build_library ALL)
#add_dependencies(build_library numerics_kernel_lib)
