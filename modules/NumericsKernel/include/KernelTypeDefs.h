#include <blaze/Blaze.h>

namespace KERNEL {
    using vector = blaze::DynamicVector<scalar, blaze::columnVector>;
    using smatrix = blaze::CompressedMatrix<scalar, blaze::rowMajor>;
    using dmatrix = blaze::DynamicMatrix<scalar, blaze::rowMajor>;
}
