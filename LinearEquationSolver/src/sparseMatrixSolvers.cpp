#include "sparseMatrixSolvers.h"
#include <iostream>

// to do, move this function to another place
blaze::DynamicMatrix<double> inline toBlaze(const std::vector<std::vector<double>>& v) {
    if( v.empty() ) return blaze::DynamicMatrix<double>();   // tom matrix

    const size_t m = v.size();
    const size_t n = v[0].size();

    blaze::DynamicMatrix<double> A(m, n);

    for (size_t i = 0; i < m; ++i) {

        for (size_t j = 0; j < n; ++j) {
            A(i,j) = v[i][j];
        }
    }

    return A;
}

// to do, move this function to another place
blaze::DynamicVector<double> inline toBlaze(const std::vector<double>& v) {
    blaze::DynamicVector<double> b( v.size() );
    std::copy(v.begin(), v.end(), b.begin());
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

    const double norm_b = std::max(norm(b_), 1e-30);
    double norm_res = norm(r0);
    double rel_res  = norm_res / norm_b;

    std::size_t it = 0;
    while (rel_res > tolerance_ && it < maxIter_)
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

        norm_res = norm(r);
        rel_res  = norm_res / norm_b;
        ++it;
    }

    if (rel_res <= tolerance_)
        std::cout << "Bi_CG_STAB_sparse converged in " << it << " iterations.\n";
    else
        std::cout << "Bi_CG_STAB_sparse NOT converged (iters=" << it << ", rel res=" << rel_res << ").\n";


    return {x.begin(), x.end()};
}