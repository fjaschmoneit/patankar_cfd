#include <vector>
#include <blaze/Math.h>
#include "TypeDefs_NumericsKernel.h"

namespace Sparse{

    class matrix{

    public:

        blaze::CompressedMatrix<LINALG::scalar> data_;
        // blaze::DynamicVector<LINALG::scalar> b_;
        // blaze::DynamicVector<LINALG::scalar> x0_;
        blaze::Band<LINALG::matrix> ap_;
        blaze::Band<LINALG::matrix> ae_;
        blaze::Band<LINALG::matrix> aw_;
        blaze::Band<LINALG::matrix> an_;
        blaze::Band<LINALG::matrix> as_;

        // LINALG::scalar tolerance_ = 1e-15;
        // LINALG::scalar maxIter_ = 10000000;
        std::vector<int> bands_;

        // LINALG::vector b_;
        // int meshSize_; // why?

        // ISparseLinearSolver(const std::vector<std::vector<LINALG::scalar>> &A, const LINALG::vector &b,int meshSize);

        matrix(unsigned int nx, unsigned int ny);

        // ISparseLinearSolver(int meshSize);


        // virtual ~ISparseLinearSolver() = default;
        // virtual std::vector<LINALG::scalar> solve() = 0;

        // dir = 0:centre; 1:east; 2:north, 3:west, 4:south
        void setDirectionalFlux( const std::vector<LINALG::scalar>& ai,const unsigned int dir );
    };





    void solve_BiCGSTAB(
        LINALG::matrix&  A,
        LINALG::vector& b,
        LINALG::vector& x,      // initial guess comes in here, result is replacing that guess
        LINALG::scalar tolerance,
        unsigned int maxIter );

    //
    // class GaussSeidel : public ISparseLinearSolver {
    //
    // public:
    //     // using base class constructor:
    //     using ISparseLinearSolver::ISparseLinearSolver;
    //
    //     std::vector<LINALG::scalar> solve() override;
    // };
    //
    // class BiCGSTAB : public ISparseLinearSolver {
    //
    // public:
    //     // using base class constructor:
    //     using ISparseLinearSolver::ISparseLinearSolver;
    //
    //     std::vector<LINALG::scalar> solve() override;
    // };


}
