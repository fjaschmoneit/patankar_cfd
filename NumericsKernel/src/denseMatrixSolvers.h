
#include <vector>
#include <blaze/Math.h>
#include "TypeDefs_NumericsKernel.h"


namespace Dense {

    // virtual class defining the interface for all dense solvers
    class IDenseLinearSolver {
    protected:

        blaze::DynamicMatrix<LINALG::scalar> A_;
        blaze::DynamicVector<LINALG::scalar> b_;
        blaze::DynamicVector<LINALG::scalar> x0_;
        LINALG::scalar tolerance_ = 1e-14;
        LINALG::scalar maxIter_ = 100000;


    public:
        void setX0(std::vector<LINALG::scalar> x0);

        IDenseLinearSolver(const std::vector<std::vector<LINALG::scalar>> &A, const std::vector<LINALG::scalar> &b);
        //      : x0_(b_.size());

        virtual ~IDenseLinearSolver() = default;

        virtual std::vector<LINALG::scalar> solve() = 0;
    };

    class GaussSeidel : public IDenseLinearSolver
    {
        public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<LINALG::scalar>  solve() override ;
    };
    class JacobiIter : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<LINALG::scalar> solve() override ;
    };

    class BiCGSTAB : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<LINALG::scalar>  solve() override ;
    };

}
