#include "LinEqsSolvers.h"

template<typename MatrixType>
matrix<MatrixType>::matrix(const unsigned int nx, const unsigned int ny):
    data_(nx*ny, nx*ny),
    bands_(std::vector<int>{ static_cast<int>(-nx), -1, 0, 1, static_cast<int>(nx)} ),
    ap_(blaze::band(data_,0)),
    an_(blaze::band(data_,static_cast<int>(-nx))),
    ae_(blaze::band(data_,1)),
    as_(blaze::band(data_,static_cast<int>(nx))),
    aw_(blaze::band(data_,-1))
{
}



// void Dense::IDenseLinearSolver::setX0(std::vector<LINALG::scalar> x0)
// {
//      blaze::reset(x0_);
//
//      x0_ =  toBlaze(x0);
// }


// std::vector<LINALG::scalar> Dense::JacobiIter::solve() {
//
//     auto rows = A_.rows();
//     blaze::DynamicVector<LINALG::scalar> x(rows, 0.0 );
//     blaze::DynamicVector<LINALG::scalar> x_old = x0_;
//     blaze::DynamicMatrix<LINALG::scalar> B = A_;     // deep copy
//
//     auto d = blaze::diagonal(A_);  // view af diagonalen
//     auto invD = 1.0/d;
//
//     blaze::DiagonalMatrix< blaze::CompressedMatrix<LINALG::scalar> > invDSparse( A_.rows() );
//     blaze::diagonal( invDSparse ) = invD;
//     blaze::diagonal(B) = 0.0;
//
//     for (int k = 0; k < maxIter_; ++k)
//     {
//         auto rhs = b_ - B * x_old;
//         x = invDSparse * rhs;
//         if ( MathUtils::diffNorm2(x,x_old) < tolerance_)
//         {
//             std::cout<<"Jacobi solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
//             return std::vector<double> (x.begin(), x.end());
//         }
//         x_old = x;
//     }
//     std::cerr<<"Jacobi solver did not converge within"<<  std::to_string(maxIter_)<<" iterations."<<std::endl;
//     return {x.begin(), x.end()};
// }


// std::vector<LINALG::scalar> Dense::BiCGSTAB::solve()
// {
//
//     using blaze::dot;
//     using blaze::norm;
//
//     const std::size_t n = A_.rows();
//     blaze::DynamicVector<LINALG::scalar> x = x0_;
//     blaze::DynamicVector<LINALG::scalar> r0 = b_ - A_ * x;
//     blaze::DynamicVector<LINALG::scalar> r  = r0;
//     blaze::DynamicVector<LINALG::scalar> p  = r;
//     blaze::DynamicVector<LINALG::scalar> v(n, 0.0), s(n, 0.0), t(n, 0.0);
//
//     LINALG::scalar rho  = dot(r0, r0);
//     LINALG::scalar alpha = 0.0, omega = 0.0, rho1 = 0.0, beta = 0.0;
//
//     const LINALG::scalar normb = std::max(norm(b_), 1e-30);
//     LINALG::scalar normres = norm(r0);
//     LINALG::scalar relres  = normres / normb;
//
//     std::size_t it = 0;
//     while (relres > tolerance_ && it < maxIter_)
//     {
//         v    = A_ * p;
//         LINALG::scalar vr0 = dot(v, r0);
//         if (std::fabs(vr0) < 1e-30) break;
//
//         alpha = rho / vr0;
//
//         s = r - alpha * v;
//         t = A_ * s;
//
//         LINALG::scalar tt = dot(t, t);
//         if (tt <= 0.0) break;
//
//         omega = dot(t, s) / tt;
//         if (std::fabs(omega) < 1e-30) break;
//
//         x = x + alpha * p + omega * s;
//         r = s - omega * t;
//
//         rho1 = dot(r, r0);
//         beta = (rho1 / rho) * (alpha / omega);
//         p = r + beta * (p - omega * v);
//         rho = rho1;
//
//         normres = norm(r);
//         relres  = normres / normb;
//         ++it;
//     }
//     if (relres <= tolerance_)
//         std::cout << "BiCGSTAB converged in " << it << " iterations.\n";
//     else
//         std::cout << "BiCGSTAB NOT converged (iters=" << it << ", relres=" << relres << ").\n";
//
//
//     return {x.begin(), x.end()};
// }
