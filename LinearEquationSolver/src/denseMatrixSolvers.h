
#include <vector>
#include <blaze/Math.h>

namespace Dense {

    // virtual class defining the interface for all dense solvers
    class IDenseLinearSolver {

    protected:

        blaze::DynamicMatrix<double> Ab_;
        // to do, temporary function, will be removed later
        std::vector<std::vector<double>> A_;
        // to do, temporary function, will be removed later
        std::vector<double> b_;
        blaze::DynamicVector<double> bb_;
        double tolerance_ = 1e-18;
        double maxIter_ = 10000000;


    public:
        // to do, temporary function, will be removed later
        IDenseLinearSolver(const std::vector<std::vector<double>> &A, const std::vector<double> &b)
        : A_(A), b_(b) {}
        // to do, temporary function, will be removed later
        IDenseLinearSolver(const blaze::DynamicMatrix<double> &A, const blaze::DynamicVector<double> &b)
        : Ab_(A), bb_(b) {}

        virtual ~IDenseLinearSolver() = default;

        virtual std::vector<double> solve() = 0;
        virtual blaze::DynamicVector<double> solve1() = 0;
    };

    class GaussSeidel : public IDenseLinearSolver
    {
        public:
        using IDenseLinearSolver::IDenseLinearSolver;
        // to do, temporary function, will be removed later
        std::vector<double> solve() override ;
        blaze::DynamicVector<double> solve1() override ;
    };
    class JacobiIter : public IDenseLinearSolver
    {
    public:
        using IDenseLinearSolver::IDenseLinearSolver;
        // to do, temporary function, will be removed later
        std::vector<double> solve() override ;
        blaze::DynamicVector<double> solve1() override ;
    };

}
