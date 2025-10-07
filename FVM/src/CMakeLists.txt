
# FJA: I am uncertain what this is.
# Create static library for the linear solver
add_library(modelEquation_lib STATIC
            FVMCalculateCoefficient.cpp
            FVMCalculateCoefficient.h
            )

# Set include directories for the library
target_include_directories(modelEquation_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

include_directories(/usr/local/include)
