
add_library(linear_solver_lib STATIC IMPORTED)
set_target_properties(linear_solver_lib PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/../NumericsKernel/cmake-build-debug/src/liblinear_solver_lib.a"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../NumericsKernel/src/"
)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

