# FVM/CMakeLists.txt


## Choose precision here
## For float:
#target_compile_definitions(linear_solver_lib PUBLIC KERNEL_SCALAR_T=float)
## Or for double:
## target_compile_definitions(kernel_lib PUBLIC KERNEL_SCALAR_T=double)

# Add subdirectories
add_subdirectory(src)

