#include "sparseMatrixSolvers.h"
#include <iostream>
#include <algorithm>

//
// Sparse::ISparseLinearSolver::ISparseLinearSolver(const std::vector<std::vector<LINALG::scalar>> &A, const LINALG::vector &b,int meshSize):
//      A_(toBlaze(A)),
//      b_(b),
//      bands_(std::vector<int>{ -meshSize-1, -1, 0, 1, meshSize+1} ),
//      ap_(blaze::band(A_,0)),
//      an_(blaze::band(A_,-meshSize)),
//      ae_(blaze::band(A_,1)),
//      as_(blaze::band(A_,meshSize)),
//      aw_(blaze::band(A_,-1))
// {
//     x0_.resize(b.size());
// }


Sparse::matrix::matrix(const unsigned int nx, const unsigned int ny):
    data_(nx*ny, nx*ny),
    bands_(std::vector<int>{ static_cast<int>(-nx), -1, 0, 1, static_cast<int>(nx)} ),
    ap_(blaze::band(data_,0)),
    an_(blaze::band(data_,static_cast<int>(-nx))),
    ae_(blaze::band(data_,1)),
    as_(blaze::band(data_,static_cast<int>(nx))),
    aw_(blaze::band(data_,-1))
{
    // x0_.resize(nx*ny);
}

//
// Sparse::ISparseLinearSolver::ISparseLinearSolver(unsigned int nx, unsigned int ny):
//     A_( ny*nx, ny*nx ),
//     bands_(std::vector<int>{ static_cast<int>(-nx), -1, 0, 1, static_cast<int>(nx)} ),
//     ap_(blaze::band(A_,0)),
//     an_(blaze::band(A_,static_cast<int>(-nx))),
//     ae_(blaze::band(A_,1)),
//     as_(blaze::band(A_,static_cast<int>(nx))),
//     aw_(blaze::band(A_,-1)),
//     b_(nx*ny,0),
//     meshSize_(nx*ny)
// {
//     x0_.resize(meshSize_);
// }



//
// // I don;t get it...
// Sparse::ISparseLinearSolver::ISparseLinearSolver(int meshSize):
//      A_( meshSize*meshSize, meshSize*meshSize ),
//      bands_(std::vector<int>{ static_cast<int>(-meshSize), -1, 0, 1, static_cast<int>(meshSize)} ),
//      ap_(blaze::band(A_,0)),
//      an_(blaze::band(A_,static_cast<int>(-meshSize))),
//      ae_(blaze::band(A_,1)),
//      as_(blaze::band(A_,static_cast<int>(meshSize))),
//      aw_(blaze::band(A_,-1)),
//      b_(meshSize*meshSize,0),
//     meshSize_(meshSize)
// {
//     x0_.resize(meshSize*meshSize);
// }

//
// std::vector<LINALG::scalar> Sparse::GaussSeidel::solve()
// {
//     std::vector<LINALG::scalar> x;
//     if (A_.rows() == 0)
//     {
//         x.resize(0);
//         return x;
//     }
//     return {x.begin(), x.end()};
// }

void Sparse::matrix::setDirectionalFlux( const std::vector<LINALG::scalar>& ai, const unsigned int dir )
{
    auto copyElements = [](blaze::Band<LINALG::matrix> & matrixBand, std::vector<LINALG::scalar> const & fullVec, size_t first)
    {
        for(auto i=0; i<matrixBand.size(); ++i)
        {
            matrixBand[i] = fullVec[ i + first ];
        }
    };

    if (dir == 0) copyElements( ap_, ai, 0 );
    else if (dir == 1) copyElements( ae_, ai, 0 );
    else if (dir == 2) copyElements( an_, ai, 0 );
    else if (dir == 3) copyElements( aw_, ai, 0 );
    else if (dir == 4) copyElements( as_, ai, 0 );

}

void Sparse::solve_BiCGSTAB(LINALG::matrix &A, LINALG::vector& b,
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


    const LINALG::scalar normb = std::max(norm( b ), 1e-30);  // ops!!! I shouldn't use this copy too often
    LINALG::scalar normres = norm(r0);
    LINALG::scalar relres  = normres / normb;

    const double norm_b = std::max(norm( b ), 1e-30);

    double norm_res = norm(r0);
    double rel_res  = norm_res / norm_b;

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