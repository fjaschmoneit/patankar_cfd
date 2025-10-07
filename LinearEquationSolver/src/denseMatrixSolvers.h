
#include <vector>
#include <blaze/Math.h>
#include "globalTypeDefs.h"


namespace Dense {

    // virtual class defining the interface for all dense solvers
    class IDenseLinearSolver {
    protected:

        blaze::DynamicMatrix<GLOBAL::scalar> A_;
        blaze::DynamicVector<GLOBAL::scalar> b_;
        blaze::DynamicVector<GLOBAL::scalar> x0_;
        GLOBAL::scalar tolerance_ = 1e-14;
        GLOBAL::scalar maxIter_ = 100000;


    public:
        void setX0(std::vector<GLOBAL::scalar> x0);

        IDenseLinearSolver(const std::vector<std::vector<GLOBAL::scalar>> &A, const std::vector<GLOBAL::scalar> &b);
        //      : x0_(b_.size());

        virtual ~IDenseLinearSolver() = default;

        virtual std::vector<GLOBAL::scalar> solve() = 0;
    };

    class GaussSeidel : public IDenseLinearSolver
    {
        public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<GLOBAL::scalar>  solve() override ;
    };
    class JacobiIter : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<GLOBAL::scalar> solve() override ;
    };

    class BiCGSTAB : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<GLOBAL::scalar>  solve() override ;
    };

}
