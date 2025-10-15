//#pragma once
#include <vector>
#include <blaze/Math.h>

namespace  MathUtils
{

    template <typename VectorType>
    static double norm2(const VectorType& x)
    {
        auto dq = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
        return std::sqrt(dq);
    }

    // ||x - y||_2  (Euclidean norm of the difference)
    template <typename VectorType>
    static double diffNorm2(const VectorType& x, const VectorType& y)
    {
        VectorType diff(x.size());
        std::transform(x.begin(), x.end(), y.begin(), diff.begin(),
                       [](double a, double b){ return a - b; });
        return norm2(diff);
    }


    inline blaze::DynamicMatrix<LINALG::scalar> toBlaze(const std::vector<std::vector<LINALG::scalar>>& v) {
        if( v.empty() ) return blaze::DynamicMatrix<LINALG::scalar>();   // tom matrix

        const size_t m = v.size();
        const size_t n = v[0].size();

        blaze::DynamicMatrix<LINALG::scalar> A(m, n);

        for (size_t i = 0; i < m; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                A(i,j) = v[i][j];
            }
        }

        return A;
    }

    inline blaze::DynamicVector<LINALG::scalar> toBlaze(const std::vector<LINALG::scalar>& v) {
        blaze::DynamicVector<LINALG::scalar> b( v.size() );
        std::copy(v.begin(), v.end(), b.begin());
        return b;
    }

    inline void takeUpperOrLowerStrictlyTriangular(blaze::DynamicMatrix<LINALG::scalar>& T, bool isLowerTriangular)
    {
        //lambda med sign function and k as input
        auto sign = [](bool b, auto k) {return b ? k : -k;};

        blaze::reset( blaze::diagonal( T ) );
        const long n = static_cast<long>( T.rows() );
        for( long k = 1; k <= (n-1); ++k )
        {
            blaze::reset( blaze::band( T, sign(isLowerTriangular, k) ) );
        }

    }

    inline blaze::DynamicMatrix<LINALG::scalar> copyStrictlyUpperTriangular(const blaze::DynamicMatrix<LINALG::scalar>& A)
    {
        blaze::DynamicMatrix<LINALG::scalar> R = A;
        takeUpperOrLowerStrictlyTriangular(R ,false);

        return R;
    }

    inline blaze::DynamicMatrix<LINALG::scalar> copyStrictlyLowerTriangular(const blaze::DynamicMatrix<LINALG::scalar>& A) {
        blaze::DynamicMatrix<LINALG::scalar> R = A;
        takeUpperOrLowerStrictlyTriangular(R ,true);

        return R;
    }
};


