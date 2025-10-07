#include <vector>
#include <blaze/Math.h>
#include "globalTypeDefs.h"
namespace Sparse{

    // interface
    class ISparseLinearSolver {

    protected:

        blaze::CompressedMatrix<GLOBAL::scalar> A_;
//        blaze::DynamicVector<GLOBAL::scalar> b_;
        blaze::DynamicVector<GLOBAL::scalar> x0_;
        GLOBAL::scalar tolerance_ = 1e-15;
        GLOBAL::scalar maxIter_ = 10000000;
        std::vector<int> bands_;

        blaze::Band<LINALG::matrix> ap_;
        blaze::Band<LINALG::matrix> ae_;
        blaze::Band<LINALG::matrix> aw_;
        blaze::Band<LINALG::matrix> an_;
        blaze::Band<LINALG::matrix> as_;
        //LINALG::vector b_;
        GLOBAL::vector b_;
        int meshSize_;


    public:
        ISparseLinearSolver(const std::vector<std::vector<GLOBAL::scalar>> &A, const GLOBAL::vector &b,int meshSize);

        explicit ISparseLinearSolver(int meshSize);

        virtual ~ISparseLinearSolver() = default;

        virtual std::vector<GLOBAL::scalar> solve() = 0;

        void setDirectionalFlux( const GLOBAL::vector& ai,const FVM::CardinalDirection dir );
    };

    class GaussSeidel : public ISparseLinearSolver {

    public:
        // using base class constructor:
        using ISparseLinearSolver::ISparseLinearSolver;

        std::vector<GLOBAL::scalar> solve() override;
    };

    class BiCGSTAB : public ISparseLinearSolver {

    public:
        // using base class constructor:
        using ISparseLinearSolver::ISparseLinearSolver;

        std::vector<GLOBAL::scalar> solve() override;
    };


}
