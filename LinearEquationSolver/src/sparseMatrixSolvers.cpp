#include "sparseMatrixSolvers.h"
#include <iostream>

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

Sparse::ISparseLinearSolver::ISparseLinearSolver(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
{
    A_ = toBlaze(A);
    b_ = toBlaze(b);
    x0_.resize(b.size());

}

std::vector<double> Sparse::GaussSeidel::solve()
{
    std::vector<double> x;
    if (A_.rows() == 0)
    {
        x.resize(0);
        return x;
    }



    return {x.begin(), x.end()};
}

std::vector<double> Sparse::BiCGSTAB::solve()
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
        std::cout << "BiCGSTAB_sparse converged in " << it << " iterations.\n";
    else
        std::cout << "BiCGSTAB_sparse NOT converged (iters=" << it << ", relres=" << relres << ").\n";


    return {x.begin(), x.end()};
}