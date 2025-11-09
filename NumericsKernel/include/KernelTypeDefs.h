#include <blaze/Blaze.h>

namespace KERNEL {
#ifndef KERNEL_SCALAR_T
    using scalar = double;              // default
#else
    using scalar = KERNEL_SCALAR_T;     // define this in cmake file
#endif

    using vector = blaze::DynamicVector<scalar, blaze::columnVector>;
    using smatrix = blaze::CompressedMatrix<scalar, blaze::rowMajor>;
    using dmatrix = blaze::DynamicMatrix<scalar, blaze::rowMajor>;
}
