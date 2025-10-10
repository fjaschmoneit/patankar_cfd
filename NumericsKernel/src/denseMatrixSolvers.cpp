
#include "denseMatrixSolvers.h"

#include <iostream>
#include "MathUtils.h"

blaze::DynamicMatrix<LINALG::scalar> toBlaze(const std::vector<std::vector<LINALG::scalar>>& v) {
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

  blaze::DynamicVector<LINALG::scalar> toBlaze(const std::vector<LINALG::scalar>& v) {
      blaze::DynamicVector<LINALG::scalar> b( v.size() );
      std::copy(v.begin(), v.end(), b.begin());
      return b;
  }


Dense::IDenseLinearSolver::IDenseLinearSolver(const std::vector<std::vector<LINALG::scalar> > &A,
                                              const std::vector<LINALG::scalar> &b) {
    A_ = toBlaze(A);
    b_ = toBlaze(b);
    x0_.resize(b.size());
}


void takeUpperOrLowerStrictlyTriangular(blaze::DynamicMatrix<LINALG::scalar>& T, bool isLowerTriangular)
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




blaze::DynamicMatrix<LINALG::scalar> copyStrictlyUpperTriangular(const blaze::DynamicMatrix<LINALG::scalar>& A)
{
    blaze::DynamicMatrix<LINALG::scalar> R = A;
    takeUpperOrLowerStrictlyTriangular(R ,false);

    return R;
}

blaze::DynamicMatrix<LINALG::scalar> copyStrictlyLowerTriangular(const blaze::DynamicMatrix<LINALG::scalar>& A)
{
    blaze::DynamicMatrix<LINALG::scalar> R = A;
    takeUpperOrLowerStrictlyTriangular(R ,true);

    return R;
}


std::vector<LINALG::scalar> Dense::GaussSeidel::solve()
{

    auto rows = A_.rows();
    blaze::DynamicVector<LINALG::scalar> x(rows, 0.0 );
    blaze::DynamicVector<LINALG::scalar> x_old = x0_;
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

void Dense::IDenseLinearSolver::setX0(std::vector<LINALG::scalar> x0)
{
     blaze::reset(x0_);

     x0_ =  toBlaze(x0);
}


std::vector<LINALG::scalar> Dense::JacobiIter::solve() {

    auto rows = A_.rows();
    blaze::DynamicVector<LINALG::scalar> x(rows, 0.0 );
    blaze::DynamicVector<LINALG::scalar> x_old = x0_;
    blaze::DynamicMatrix<LINALG::scalar> B = A_;     // deep copy

    auto d = blaze::diagonal(A_);  // view af diagonalen
    auto invD = 1.0/d;

    blaze::DiagonalMatrix< blaze::CompressedMatrix<LINALG::scalar> > invDSparse( A_.rows() );
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

std::vector<LINALG::scalar> Dense::BiCGSTAB::solve()
{

    using blaze::dot;
    using blaze::norm;

    const std::size_t n = A_.rows();
    blaze::DynamicVector<LINALG::scalar> x = x0_;
    blaze::DynamicVector<LINALG::scalar> r0 = b_ - A_ * x;
    blaze::DynamicVector<LINALG::scalar> r  = r0;
    blaze::DynamicVector<LINALG::scalar> p  = r;
    blaze::DynamicVector<LINALG::scalar> v(n, 0.0), s(n, 0.0), t(n, 0.0);

    LINALG::scalar rho  = dot(r0, r0);
    LINALG::scalar alpha = 0.0, omega = 0.0, rho1 = 0.0, beta = 0.0;

    const LINALG::scalar normb = std::max(norm(b_), 1e-30);
    LINALG::scalar normres = norm(r0);
    LINALG::scalar relres  = normres / normb;

    std::size_t it = 0;
    while (relres > tolerance_ && it < maxIter_)
    {
        v    = A_ * p;
        LINALG::scalar vr0 = dot(v, r0);
        if (std::fabs(vr0) < 1e-30) break;

        alpha = rho / vr0;

        s = r - alpha * v;
        t = A_ * s;

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
