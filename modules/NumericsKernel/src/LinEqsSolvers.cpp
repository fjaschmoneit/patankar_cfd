#include "LinEqsSolvers.h"
#include <blaze/math/lapack/clapack/geev.h>


namespace LINEQSOLVERS {

    template<typename MatrixType>
    GLOBAL::scalar CalMaxSpectralRadius(const MatrixType &M)
    {
        auto rows = M.rows();
        blaze::DynamicVector< std::complex<GLOBAL::scalar> > lambda( rows );
        KERNEL::dmatrix G = M; // deep copy to convert to dense matrix, as geev is not working with sparse matrices.
        blaze::geev( G, lambda );
        GLOBAL::scalar rho = blaze::max( blaze::map( lambda, [](auto c){ return std::abs(c); } ) );

        return rho;
    }

    bool solve_Jacobi(const KERNEL::dmatrix &A, KERNEL::vector& x, const KERNEL::vector &b, const GLOBAL::scalar tolerance, const unsigned int maxIter)
    {
        auto rows = A.rows();
        KERNEL::vector x_old = x;
        // We perform a deep copy to ensure that the input matrix is not altered by this function.
        // A more efficient strategy can be considered later.
        // since the A matrix is const, which is a problem if we want to altered the A matrix
        // Deep copy
        KERNEL::dmatrix B = A; // B is L+U

        //deep copy and diagonal vector as it not working for dense matrix
        KERNEL::vector dvec( A.rows() );
        dvec = blaze::diagonal( A );
        auto invD = GLOBAL::scalar(1.0) / dvec;

        blaze::DiagonalMatrix< KERNEL::smatrix > invDSparse( A.rows() );
        blaze::diagonal( invDSparse ) = invD;
        blaze::diagonal(B) = 0.0;

        // todo calc of Spectral radius is computationally heavy, as geev is not working on sparse matrix and is slowly,
        // Gershgorin circle theorem should be used instead of. See  https://en.wikipedia.org/wiki/Gershgorin_circle_theorem
        // Mj = invDSparse * (-B)  Iteration matrix
        // Nj = invDSparse * b
        // x = Mj * x_old + Nj;
        //auto rho = CalMaxSpectralRadius(invDSparse * (-B));
        int k;
        for (k = 0; k < maxIter; ++k)
        {
            auto rhs = b - B * x_old;
            // Splitting method equation x = Mj * x_old + Nj; but we divide it into
            // x = D^-1 * (b - B * X_old) to save memory, same executions time.
            x = invDSparse * rhs;
            if ( blaze::norm(x-x_old) < tolerance)
            {
                std::cout<<"Jacobi solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
                break;
            }
            x_old = x;
        }
        if ( k >= maxIter)
        {
            std::cerr<<"Jacobi solver did not converge within "<<  std::to_string(maxIter)<<" iterations."<<std::endl;
            return false;
        }
        return true;
    }

    bool solve_GaussSeidel(const KERNEL::dmatrix& A, KERNEL::vector& x, const KERNEL::vector& b, const GLOBAL::scalar tolerance, const unsigned int maxIter){

        auto n = A.rows();
        KERNEL::vector x_old = x;

        // copy
        auto DL= A;

        // DL = lower triangular with diagonal
        for(size_t i = 0; i < n; ++i) {
            for(size_t j = i+1; j < n; ++j) {
                DL(i, j) = 0.0;
            }
        }

        auto invDL = blaze::inv( DL );

        auto c = invDL*b;
        auto G = -invDL*(A-DL);

        for( int k = 0; k < maxIter; ++k )
        {
            // Mj = -invDL*(A-DL)  Iteration matrix
            // Nj =  invDL*b
            // x = Mj * x_old + Nj;
            // Using G and c for faster code.
            // x = G * x_old + c
            // x = invDL* (bb_-R * x_old);
            x = G*x_old + c;

            if( blaze::norm( x - x_old ) < tolerance )
            {
                std::cout << "Gauss-Seidel solver converged after " << k << " iterations." << std::endl;
                return true;
            }
            x_old = x;
        }
        return false;
    }
}