#include "sparseMatrixSolvers.h"
#include <iostream>
#include <algorithm>

// to do, move this function to another place
blaze::DynamicMatrix<GLOBAL::scalar> inline toBlaze(const std::vector<std::vector<GLOBAL::scalar>>& v) {
    if( v.empty() ) return blaze::DynamicMatrix<GLOBAL::scalar>();   // tom matrix



    const size_t m = v.size();
    const size_t n = v[0].size();



    blaze::DynamicMatrix<GLOBAL::scalar> A(m, n);

    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            A(i,j) = v[i][j];
        }
    }

    return A;
}

// to do, move this function to another place
blaze::DynamicVector<GLOBAL::scalar> inline toBlaze(const std::vector<GLOBAL::scalar>& v) {
    blaze::DynamicVector<GLOBAL::scalar> b( v.size() );


    std::copy(v.begin(), v.end(), b.begin());
    return b;
}

Sparse::ISparseLinearSolver::ISparseLinearSolver(const std::vector<std::vector<GLOBAL::scalar>> &A, const GLOBAL::vector &b,int meshSize):

     A_(toBlaze(A)),
     b_(b),
     bands_(std::vector<int>{ -meshSize-1, -1, 0, 1, meshSize+1} ),
     ap_(blaze::band(A_,0)),
     an_(blaze::band(A_,-meshSize)),
     ae_(blaze::band(A_,1)),
     as_(blaze::band(A_,meshSize)),
     aw_(blaze::band(A_,-1)),
     meshSize_(meshSize)
{
    x0_.resize(b.size());
}

Sparse::ISparseLinearSolver::ISparseLinearSolver(int meshSize):
     A_( meshSize*meshSize, meshSize*meshSize, 0.0 ),
     bands_(std::vector<int>{ static_cast<int>(-meshSize), -1, 0, 1, static_cast<int>(meshSize)} ),
     ap_(blaze::band(A_,0)),
     an_(blaze::band(A_,static_cast<int>(-meshSize))),
     ae_(blaze::band(A_,1)),
     as_(blaze::band(A_,static_cast<int>(meshSize))),
     aw_(blaze::band(A_,-1)),
     b_(meshSize*meshSize,0),
    meshSize_(meshSize)
{
    x0_.resize(meshSize*meshSize);
}

std::vector<GLOBAL::scalar> Sparse::GaussSeidel::solve()
{
    std::vector<GLOBAL::scalar> x;
    if (A_.rows() == 0)
    {
        x.resize(0);
        return x;
    }
    return {x.begin(), x.end()};
}

void Sparse::ISparseLinearSolver::setDirectionalFlux( const GLOBAL::vector& ai,
                                                      const FVM::CardinalDirection dir )
{
    auto copyElements = [](blaze::Band<LINALG::matrix> & matrixBand, GLOBAL::vector const & fullVec, size_t first)
    {
        for(auto i=0; i<matrixBand.size(); ++i)
        {
            matrixBand[i] = fullVec[ i + first ];
        }
    };

    switch (dir) {
        case FVM::CardinalDirection::centre:
        {
            copyElements( ap_, ai, 0 );
            break;
        }
        case FVM::CardinalDirection::east:
        {
            copyElements( ae_, ai, 0 );
			break;
        }
        case FVM::CardinalDirection::south:
        {
            copyElements( as_, ai, 0 );
			break;
        }
        case FVM::CardinalDirection::west:
        {
            copyElements( aw_, ai, 1 );
			break;
        }
        case FVM::CardinalDirection::north:
        {
            copyElements( an_, ai, meshSize_ );
			break;
        }
        case FVM::CardinalDirection::su:
        {
            b_ = GLOBAL::vector(ai.begin(), ai.end());
            break;
        }
        default: break;
    }
}

std::vector<GLOBAL::scalar> Sparse::BiCGSTAB::solve()
{
    if( blaze::isZero( A_ ) )
    {
        std::cerr << "Error in A is empty, try collect the A matrix before solve" << std::endl;
        std::vector<GLOBAL::scalar> empty;
        return {};
    }
    if (b_.size() != A_.rows())
    {
        std::cerr << "Error in B is not same side as Row of A" << std::endl;
    }
    using blaze::dot;
    using blaze::norm;

    const std::size_t n = A_.rows();
    blaze::DynamicVector<GLOBAL::scalar> x = x0_;
    blaze::DynamicVector<GLOBAL::scalar> r0 = b_.toBlaze() - A_ * x;
    blaze::DynamicVector<GLOBAL::scalar> r  = r0;
    blaze::DynamicVector<GLOBAL::scalar> p  = r;
    blaze::DynamicVector<GLOBAL::scalar> v(n, 0.0), s(n, 0.0), t(n, 0.0);

    GLOBAL::scalar rho  = dot(r0, r0);
    GLOBAL::scalar alpha = 0.0, omega = 0.0, rho1 = 0.0, beta = 0.0;


    const GLOBAL::scalar normb = std::max(norm(b_.toBlaze()), 1e-30);
    GLOBAL::scalar normres = norm(r0);
    GLOBAL::scalar relres  = normres / normb;

    const double norm_b = std::max(norm(b_.toBlaze()), 1e-30);

    double norm_res = norm(r0);
    double rel_res  = norm_res / norm_b;

    std::size_t it = 0;
    while (rel_res > tolerance_ && it < maxIter_)
    {
        v    = A_ * p;
        GLOBAL::scalar vr0 = dot(v, r0);
        if (std::fabs(vr0) < 1e-30) break;

        alpha = rho / vr0;

        s = r - alpha * v;
        t = A_ * s;

        GLOBAL::scalar tt = dot(t, t);
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

    if (rel_res <= tolerance_)
        std::cout << "Bi_CG_STAB_sparse converged in " << it << " iterations\n";
    else
        std::cout << "Bi_CG_STAB_sparse NOT converged (iters=" << it << ", rel res=" << rel_res << ").\n";



    return {x.begin(), x.end()};
}