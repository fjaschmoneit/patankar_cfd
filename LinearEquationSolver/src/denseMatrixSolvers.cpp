
#include "denseMatrixSolvers.h"

#include <iostream>
#include "MathUtils.h"


// to do, move this function to another place
blaze::DynamicMatrix<double> inline toBlaze(const std::vector<std::vector<double>>& v) {
    if( v.empty() ) return blaze::DynamicMatrix<double>();   // tom matrix

    const size_t m = v.size();       // antal rækker
    const size_t n = v[0].size();    // antal kolonner

    blaze::DynamicMatrix<double> A(m, n);

    for (size_t i = 0; i < m; ++i) {
        // evt. tjek at v[i].size() == n, ellers er input ikke rektangulært
        for (size_t j = 0; j < n; ++j) {
            A(i,j) = v[i][j];
        }
    }

    return A;
}

// to do, move this function to another place
blaze::DynamicVector<double> inline toBlaze(const std::vector<double>& v) {
    blaze::DynamicVector<double> b( v.size() );           // allokér
    std::copy(v.begin(), v.end(), b.begin());             // kopiér
    return b;
}

Dense::IDenseLinearSolver::IDenseLinearSolver(const std::vector<std::vector<double> > &A,
                                              const std::vector<double> &b) {
    A_ = toBlaze(A);
    b_ = toBlaze(b);
    x0_.resize(b.size());
}


void takeUpperOrLowerStrictlyTriangular(blaze::DynamicMatrix<double>& T, bool isLowerTriangular)
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



blaze::DynamicMatrix<double> copyStrictlyUpperTriangular(const blaze::DynamicMatrix<double>& A)
{
    blaze::DynamicMatrix<double> R = A;
    takeUpperOrLowerStrictlyTriangular(R ,false);

    return R;
}

blaze::DynamicMatrix<double> copyStrictlyLowerTriangular(const blaze::DynamicMatrix<double>& A)
{
    blaze::DynamicMatrix<double> R = A;
    takeUpperOrLowerStrictlyTriangular(R ,true);

    return R;
}


std::vector<double> Dense::GaussSeidel::solve()
{

    auto rows = A_.rows();
    blaze::DynamicVector<double> x(rows, 0.0 );
    blaze::DynamicVector<double> x_old = x0_;
    auto R = copyStrictlyUpperTriangular(A_);
    auto L = copyStrictlyLowerTriangular(A_);
    blaze::DynamicMatrix<double> D( A_.rows(), A_.columns(), 0.0 );
    blaze::diagonal(D,0) = blaze::diagonal(A_,0);
    auto DL = D + L;  // (D+L)
    auto invDL = blaze::inv( DL );
    auto c = invDL*b_;
    auto G = -invDL*R;

    for( int k = 0; k < maxIter_; ++k )
    {
        // Using G and c for faster code.
        // x = G * x_old + c
        // x = invDL* (bb_-R * x_old);
        x = G*x_old + c;
        if( MathUtils::diffNorm2( x, x_old ) < tolerance_ )
        {
            std::cout << "Gauss-Seidel solver converged after " << k << " iterations." << std::endl;
            return std::vector<double> (x.begin(), x.end());
        }
        x_old = x;
    }

    return {x.begin(), x.end()};
}

void Dense::IDenseLinearSolver::setX0(std::vector<double> x0)
{
     blaze::reset(x0_);

     x0_ =  toBlaze(x0);
}


std::vector<double> Dense::JacobiIter::solve() {

    auto rows = A_.rows();
    blaze::DynamicVector<double> x(rows, 0.0 );
    blaze::DynamicVector<double> x_old = x0_;
    blaze::DynamicMatrix<double> B = A_;     // deep copy

    auto d = blaze::diagonal(A_);  // view af diagonalen
    auto invD = 1.0/d;

    blaze::DiagonalMatrix< blaze::CompressedMatrix<double> > invDSparse( A_.rows() );
    blaze::diagonal( invDSparse ) = invD;
    blaze::diagonal(B) = 0.0;

    for (int k = 0; k < maxIter_; ++k)
    {
        auto rhs = b_ - B * x_old;
        x = invDSparse * rhs;
        if ( MathUtils::diffNorm2(x,x_old) < tolerance_)
        {
            std::cout<<"Jacobi solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
            return std::vector<double> (x.begin(), x.end());
        }
        x_old = x;
    }
    std::cerr<<"Jacobi solver did not converge within"<<  std::to_string(maxIter_)<<" iterations."<<std::endl;
    return {x.begin(), x.end()};
}

std::vector<double> Dense::BiCGSTAB::solve()
{

    using blaze::dot;
    using blaze::norm;

    const std::size_t n = A_.rows();
    blaze::DynamicVector<double> x = x0_;
    blaze::DynamicVector<double> r0 = b_ - A_ * x;
    blaze::DynamicVector<double> r  = r0;
    blaze::DynamicVector<double> p  = r;
    blaze::DynamicVector<double> v(n, 0.0), s(n, 0.0), t(n, 0.0);

    double rho  = dot(r0, r0);
    double alpha = 0.0, omega = 0.0, rho1 = 0.0, beta = 0.0;

    const double normb = std::max(norm(b_), 1e-30);
    double normres = norm(r0);
    double relres  = normres / normb;

    std::size_t it = 0;
    while (relres > tolerance_ && it < maxIter_)
    {
        v    = A_ * p;
        double vr0 = dot(v, r0);
        if (std::fabs(vr0) < 1e-30) break;

        alpha = rho / vr0;

        s = r - alpha * v;
        t = A_ * s;

        double tt = dot(t, t);
        if (tt <= 0.0) break;

        omega = dot(t, s) / tt;
        if (std::fabs(omega) < 1e-30) break;

        x = x + alpha * p + omega * s;
        r = s - omega * t;

        rho1 = dot(r, r0);
        beta = (rho1 / rho) * (alpha / omega);
        p = r + beta * (p - omega * v);
        rho = rho1;

        normres = norm(r);
        relres  = normres / normb;
        ++it;
    }

    if (relres <= tolerance_)
        std::cout << "BiCGSTAB converged in " << it << " iterations.\n";
    else
        std::cout << "BiCGSTAB NOT converged (iters=" << it << ", relres=" << relres << ").\n";


    return {x.begin(), x.end()};
}
