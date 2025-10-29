#include <vector>
#include <blaze/Math.h>
#include "TypeDefs_NumericsKernel.h"
#include <iostream>
#include "MathUtils.h"
#include <cmath>



template<typename MatrixType>
class matrix {

public:
    MatrixType data_;
    blaze::Band<MatrixType> ap_;
    blaze::Band<MatrixType> ae_;
    blaze::Band<MatrixType> aw_;
    blaze::Band<MatrixType> an_;
    blaze::Band<MatrixType> as_;
    std::vector<int> bands_;

    matrix(unsigned int nx, unsigned int ny);

};

template<typename MatrixType>
void fillBand(blaze::Band<MatrixType> band, LINALG::scalar value) {
    for(size_t i = 0; i < band.size(); i++)  band[i] = value;
}

template<typename MatrixType>
void fillBand(blaze::Band<MatrixType> band, LINALG::vector values) {
    for(size_t i = 0; i < band.size(); i++)  band[i] = values[i];
}

//
// template<typename MatrixType>
// void solve_BiCGSTAB(
//     MatrixType&  A,
//     LINALG::vector& b,
//     LINALG::vector& x,      // initial guess comes in here, result is replacing that guess
//     LINALG::scalar tolerance,
//     unsigned int maxIter );

// free template class I cannot separate definition from implementation. That leads to linking error.
template<typename MatrixType>
void solve_BiCGSTAB(MatrixType &A, LINALG::vector& b,
    LINALG::vector& x, LINALG::scalar tolerance, unsigned int maxIter) {

    if( blaze::isZero( A ) )
    {
        throw std::invalid_argument("Error in A is empty");
    }
    if (b.size() != A.rows())
    {
        throw std::invalid_argument("Error in B is not same side as Row of A");
    }
    using blaze::dot;
    using blaze::norm;

    const std::size_t n = A.rows();
    LINALG::vector r0( b - A * x );
    LINALG::vector r  = r0;
    LINALG::vector p  = r;
    LINALG::vector v(n, 0.0), s(n, 0.0), t(n, 0.0);

    LINALG::scalar rho  = dot(r0, r0);
    LINALG::scalar alpha = 0.0, omega = 0.0, rho1 = 0.0, beta = 0.0;


    const LINALG::scalar normb = std::max(norm( b ), 1e-30);
    LINALG::scalar normres = norm(r0);
    LINALG::scalar relres  = normres / normb;

    const LINALG::scalar norm_b = std::max(norm( b ), 1e-30);

    LINALG::scalar norm_res = norm(r0);
    LINALG::scalar rel_res  = norm_res / norm_b;

    std::size_t it = 0;
    while (rel_res > tolerance && it < maxIter)
    {
        v    = A * p;
        LINALG::scalar vr0 = dot(v, r0);
        if (std::fabs(vr0) < 1e-30) break;

        alpha = rho / vr0;

        s = r - alpha * v;
        t = A * s;

        LINALG::scalar tt = dot(t, t);
        if (tt <= 0.0) break;

        omega = dot(t, s) / tt;
        if (std::fabs(omega) < 1e-30) break;

        x = x + alpha * p + omega * s;
        r = s - omega * t;

        rho1 = dot(r, r0);
        beta = (rho1 / rho) * (alpha / omega);
        p = r + beta * (p - omega * v);
        rho = rho1;

        norm_res = norm(r);
        rel_res  = norm_res / norm_b;
        ++it;
    }

    if (rel_res <= tolerance)
        std::cout << "Bi_CG_STAB_sparse converged in " << it << " iterations\n";
    else
        std::cout << "Bi_CG_STAB_sparse NOT converged (iters=" << it << ", rel res=" << rel_res << ").\n";
}

template<typename MatrixType>
void solve_GaussSeidel(MatrixType &A, LINALG::vector &b, LINALG::vector &x, LINALG::scalar tolerance, unsigned int maxIter) {
    auto rows = A.rows();
    // blaze::DynamicVector<LINALG::scalar> x(rows, 0.0 );
    blaze::DynamicVector<LINALG::scalar> x_old = x;
    auto R = MathUtils::copyStrictlyUpperTriangular(A);            // can we not do better with a native balze function?
    auto L = MathUtils::copyStrictlyLowerTriangular(A);
    blaze::DynamicMatrix<double> D( A.rows(), A.columns(), 0.0 );
    blaze::diagonal(D,0) = blaze::diagonal(A,0);
    auto DL = D + L;  // (D+L)
    auto invDL = blaze::inv( DL );
    auto c = invDL*b;
    auto G = -invDL*R;

    for( int k = 0; k < maxIter; ++k )
    {
        // Using G and c for faster code.
        // x = G * x_old + c
        // x = invDL* (bb_-R * x_old);
        x = G*x_old + c;

        // something's wrong here
        if( MathUtils::diffNorm2( x, x_old ) < tolerance )
        {
            std::cout << "Gauss-Seidel solver converged after " << k << " iterations." << std::endl;
            break;
        }
        x_old = x;
    }
}

template<typename MatrixType>
void solve_JacobiIter(MatrixType &A, LINALG::vector &b, LINALG::vector &x, LINALG::scalar tolerance, unsigned int maxIter)
{

     auto rows = A.rows();
     //blaze::DynamicVector<LINALG::scalar> x_old ;
     // we need to make a xold as input.
     blaze::DynamicVector<double> x_old( rows, 0.0 );
     blaze::DynamicMatrix<LINALG::scalar> B = A;     // deep copy



    //deep copy and diagonal vector as it not working for dense matrix
     blaze::DynamicVector<LINALG::scalar> dvec( A.rows() );
     dvec = blaze::diagonal( A );
     auto invD = LINALG::scalar(1.0) / dvec;

     //for dense Matrix, used the this
     //auto d = blaze::diagonal(A);  // view af diagonalen
     //auto invD = 1.0/d;

     blaze::DiagonalMatrix< blaze::CompressedMatrix<LINALG::scalar> > invDSparse( A.rows() );
     blaze::diagonal( invDSparse ) = invD;
     blaze::diagonal(B) = 0.0;
    int k;
     for (k = 0; k < maxIter; ++k)
     {
         auto rhs = b - B * x_old;
         x = invDSparse * rhs;
         if ( MathUtils::diffNorm2(x,x_old) < tolerance)
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





//
// void solve_Jacobi(
//     LINALG::matrix&  A,
//     LINALG::vector& b,
//     LINALG::vector& x,
//     LINALG::scalar tolerance,
//     unsigned int maxIter );
//
// void solve_GaussSeidel(
//     LINALG::matrix&  A,
//     LINALG::vector& b,
//     LINALG::vector& x,
//     LINALG::scalar tolerance,
//     unsigned int maxIter );

