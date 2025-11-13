#include "LinEqsSolvers.h"



namespace LINEQSOLVERS {

    void solve_Jacobi(const KERNEL::dmatrix &A, KERNEL::vector& x, const KERNEL::vector &b, const KERNEL::scalar tolerance, const unsigned int maxIter)
    {
        auto rows = A.rows();
        //KERNEL::vector x_old ;
        // we need to make a xold as input.
        KERNEL::vector x_old( rows, 0.0 );
        KERNEL::dmatrix B = A;     // deep copy

        //deep copy and diagonal vector as it not working for dense matrix
        KERNEL::vector dvec( A.rows() );
        dvec = blaze::diagonal( A );
        auto invD = KERNEL::scalar(1.0) / dvec;

        //for dense Matrix, used the this
        //auto d = blaze::diagonal(A);  // view af diagonalen
        //auto invD = 1.0/d;

        blaze::DiagonalMatrix< KERNEL::smatrix > invDSparse( A.rows() );
        blaze::diagonal( invDSparse ) = invD;
        blaze::diagonal(B) = 0.0;
        int k;
        for (k = 0; k < maxIter; ++k)
        {
            auto rhs = b - B * x_old;
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
            std::cerr<<"Jacobi solver did not converge within"<<  std::to_string(maxIter)<<" iterations."<<std::endl;
        }
    }

    void solve_GaussSeidel(const KERNEL::dmatrix& A, KERNEL::vector& x, const KERNEL::vector& b, const KERNEL::scalar tolerance, const unsigned int maxIter){

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

        // here I could free A,DL,invDL
        for( int k = 0; k < maxIter; ++k )
        {
            // Using G and c for faster code.
            // x = G * x_old + c
            // x = invDL* (bb_-R * x_old);
            x = G*x_old + c;

            if( blaze::norm( x - x_old ) < tolerance )
            {
                std::cout << "Gauss-Seidel solver converged after " << k << " iterations." << std::endl;
                break;
            }
            x_old = x;
        }
    }
}