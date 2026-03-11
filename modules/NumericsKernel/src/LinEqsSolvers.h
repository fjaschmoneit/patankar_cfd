#include <vector>
#include <blaze/Math.h>
#include <iostream>
#include <cmath>
#include "GlobalTypeDefs.h"
#include "KernelTypeDefs.h"     // should not be here FJA



template<typename MatrixType>
void checkLinEqSystemConsistency(const MatrixType& A, const KERNEL::vector& b) {
    if( blaze::isZero( A ) )
    {
        throw std::invalid_argument("Error in A is empty");
    }
    if (b.size() != A.rows())
    {
        throw std::invalid_argument("Error in B is not same side as Row of A");
    }
}


namespace LINEQSOLVERS {

    // A must be strictly diagonal dominant
    // Is very slow! check if it's okay
    void solve_GaussSeidel( const KERNEL::dmatrix& A, KERNEL::vector& x, const KERNEL::vector& b, const GLOBAL::scalar tolerance, const unsigned int maxIter);

    void solve_Jacobi(      const KERNEL::dmatrix &A, KERNEL::vector& x, const KERNEL::vector &b, const GLOBAL::scalar tolerance, const unsigned int maxIter);

    // todo BiCGSTAB is not working with single precisition problem, is have converge problems. We need to figure out why.
    // free template class I cannot separate definition from implementation. That leads to linking error.
    template<typename MatrixType>
    void solve_BiCGSTAB(const MatrixType &A, KERNEL::vector& x, const KERNEL::vector& b, const GLOBAL::scalar tolerance, const unsigned int maxIter) {

        const std::size_t n = A.rows();
        KERNEL::vector r0( b - A * x );
        KERNEL::vector r  = r0;
        KERNEL::vector p  = r;
        KERNEL::vector v(n, 0.0), s(n, 0.0), t(n, 0.0);

        GLOBAL::scalar rho  = blaze::dot(r0, r0);
        GLOBAL::scalar alpha = 0.0, omega = 0.0, rho1 = 0.0, beta = 0.0;

        GLOBAL::scalar rel_res  = blaze::norm(r0) ;

        std::size_t it = 0;
        while (rel_res > tolerance && it < maxIter)
        {
            v    = A * p;
            GLOBAL::scalar vr0 = blaze::dot(v, r0);
            if (std::fabs(vr0) < 1e-30)
            {
                auto message = std::string("BiCGSTAB breakdown: v·r0 ≈ 0 → cannot compute alpha (division by zero). Stopping.");
                throw std::runtime_error("Fatal error in BiCGSTAB solver: " + message);
            }

            alpha = rho / vr0;

            s = r - alpha * v;
            t = A * s;

            GLOBAL::scalar tt = blaze::dot(t, t);
            if (tt <= 0.0)
            {
                auto message = std::string("BiCGSTAB breakdown: t·t <= 0 → (tt must be positive). Stopping iterations.");
                throw std::runtime_error("Fatal error in BiCGSTAB solver: " + message);
            }

            omega = blaze::dot(t, s) / tt;
            if (std::fabs(omega) < 1e-30)
            {
                auto message ="BiCGSTAB breakdown: omega ≈ 0 (|omega|=" + std::to_string(std::fabs(omega))
                          + ") → division by omega would be unstable. Stopping iterations.\n";
                throw std::runtime_error("Fatal error in BiCGSTAB solver: " + message);
            }

            x = x + alpha * p + omega * s;
            r = s - omega * t;

            rho1 = blaze::dot(r, r0);
            beta = (rho1 / rho) * (alpha / omega);
            p = r + beta * (p - omega * v);
            rho = rho1;

            rel_res  = blaze::norm(r);
            ++it;
        }

        if (rel_res <= tolerance)
        {
            return;
        }
        auto message = std::string("Bi_CG_STAB_sparse NOT converged (iters=" + std::to_string(it) + ", rel res=" + std::to_string(rel_res) + ")");
        throw std::runtime_error("Fatal error in BiCGSTAB solver: " + message);
    }

}
