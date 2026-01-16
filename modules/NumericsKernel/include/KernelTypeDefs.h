#include <blaze/Blaze.h>
#include "GlobalTypeDefs.h"
namespace KERNEL {
    using vector = blaze::DynamicVector<GLOBAL::scalar, blaze::columnVector>;
    using smatrix = blaze::CompressedMatrix<GLOBAL::scalar, blaze::rowMajor>;
    using dmatrix = blaze::DynamicMatrix<GLOBAL::scalar, blaze::rowMajor>;
}
