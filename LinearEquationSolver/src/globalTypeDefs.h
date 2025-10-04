//
// Created by Peter Berg Ammundsen on 26/09/2025.
//

#include <blaze/Math.h>

#ifndef GLOBALTYPEDEFS_H
#define GLOBALTYPEDEFS_H

namespace GLOBAL
{
    using scalar = double;
    // Lightweight wrapper around std::vector to add .toBlaze() no-copy view
    struct vector : public std::vector<scalar>
    {
        using base = std::vector<scalar>;
        using base::base; // inherit ctors

        auto toBlaze()-> blaze::CustomVector<scalar,blaze::unaligned,blaze::unpadded,blaze::columnVector>
        {
            return { this->data(), this->size() };
        }

        // Const Blaze view (read-only)
        auto toBlaze() const-> blaze::CustomVector<const scalar,blaze::unaligned,blaze::unpadded,blaze::columnVector>
        {
            return { this->data(), this->size() };
        }
        // Ny helper: laver en GLOBAL::vector ud fra en Blaze-vector
        template<typename BlazeExpr>
        static vector toStd( const BlazeExpr& e ) {
            blaze::DynamicVector<scalar, blaze::columnVector> tmp = e;
            return vector( tmp.begin(), tmp.end() );
        }
    };
}


namespace LINALG
{

    using vector = blaze::DynamicVector<GLOBAL::scalar, blaze::columnVector>;
    using matrix = blaze::CompressedMatrix<GLOBAL::scalar, blaze::rowMajor>;

}
namespace FVM
{
    enum class CardinalDirection {centre, west, east, north, south, su};
}
#endif //GLOBALTYPEDEFS_H
