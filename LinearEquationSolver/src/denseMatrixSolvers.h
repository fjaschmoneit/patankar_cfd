
#include <vector>
#include <blaze/Math.h>


namespace Dense {

    // virtual class defining the interface for all dense solvers
    class IDenseLinearSolver {
    protected:


        blaze::DynamicMatrix<double> A_;
        blaze::DynamicVector<double> b_;
        blaze::DynamicVector<double> x0_;
        double tolerance_ = 1e-14;
        double maxIter_ = 100000;


    public:
        void setX0(std::vector<double> x0);

        IDenseLinearSolver(const std::vector<std::vector<double>> &A, const std::vector<double> &b);
        //      : x0_(b_.size());

        virtual ~IDenseLinearSolver() = default;

        virtual std::vector<double> solve() = 0;
    };

    class GaussSeidel : public IDenseLinearSolver
    {
        public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<double>  solve() override ;
    };
    class JacobiIter : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<double> solve() override ;
    };

    class BiCGSTAB : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;

        std::vector<double>  solve() override ;
    };

}
