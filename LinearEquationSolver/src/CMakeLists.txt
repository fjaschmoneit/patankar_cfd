# Create static library for the linear solver
add_library(linear_solver_lib STATIC
        sparseMatrixSolvers.cpp
        denseMatrixSolvers.cpp
        denseMatrixSolvers.h
        MathUtils.cpp
        MathUtils.h
)
include_directories(/usr/local/include)
# Set include directories for the library
target_include_directories(linear_solver_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
