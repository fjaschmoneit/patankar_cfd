# Create static library for the linear solver
add_library(linear_solver_lib STATIC
        sparseMatrixSolvers.cpp
        denseMatrixSolvers.cpp
        denseMatrixSolvers.h
        MathUtils.cpp
        MathUtils.h
)
# Set include directories for the library
target_include_directories(linear_solver_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


# ==========================
# BLAS/LAPACK via Accelerate
# (Solution A - macOS)
# --------------------------
# We enable LAPACK in Blaze and link Apple Accelerate which provides
# ==========================
# NOTE: We need to link against BLAS/LAPACK (e.g. Accelerate on macOS) because blaze::gesv calls dgesv_
include_directories(/usr/local/include)

# Find the Accelerate framework (macOS system BLAS/LAPACK)
find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)

# Tell Blaze that LAPACK is available and use the CBLAS header from Accelerate
# (ACCELERATE_NEW_LAPACK and BLAZE_BLAS_INCLUDE_FILE are typical for macOS)
target_compile_definitions(linear_solver_lib PRIVATE
    BLAZE_USE_LAPACK=1
    ACCELERATE_NEW_LAPACK=1
    BLAZE_BLAS_INCLUDE_FILE="cblas.h"
)

# Link the library with Accelerate so symbols like dgesv_ resolve
target_link_libraries(linear_solver_lib PRIVATE ${ACCELERATE_FRAMEWORK})
