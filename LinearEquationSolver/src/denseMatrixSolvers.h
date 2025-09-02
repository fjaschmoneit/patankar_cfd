
#include <vector>

namespace Dense {

    // virtual class defining the interface for all dense solvers
    class IDenseLinearSolver {

    protected:
        std::vector<std::vector<double>> A_;
        std::vector<double> b_;
        double tolerance_ = 1e-18;
        double maxIter_ = 1000;


    public:

        IDenseLinearSolver(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
        : A_(A), b_(b) {}

        virtual ~IDenseLinearSolver() = default;

        virtual std::vector<double> solve() = 0;
    };

    class Jacobi : public IDenseLinearSolver
    {
        public:
        using IDenseLinearSolver::IDenseLinearSolver;
        std::vector<double> solve() override ;
    };

    class GaussSeidel : public IDenseLinearSolver
    {
        public:
        using IDenseLinearSolver::IDenseLinearSolver;
        std::vector<double> solve() override ;
    };
    class JacobiIter : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;
        std::vector<double> solve() override ;
    };

}
