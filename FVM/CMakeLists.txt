
add_library(numerics_kernel_lib STATIC IMPORTED)

set_target_properties(numerics_kernel_lib PROPERTIES
        IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/../NumericsKernel/cmake-build-debug/src/libnumerics_kernel_lib.a"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/../NumericsKernel/include/"
)


#
## Choose precision here
## For float:
#target_compile_definitions(linear_solver_lib PUBLIC KERNEL_SCALAR_T=float)
## Or for double:
## target_compile_definitions(kernel_lib PUBLIC KERNEL_SCALAR_T=double)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

