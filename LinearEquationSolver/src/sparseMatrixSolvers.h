
#include <vector>
#include <blaze/Math.h>
namespace Sparse{

    // interface
    class ISparseLinearSolver {

    protected:

        blaze::CompressedMatrix<double> A_;
        blaze::DynamicVector<double> b_;
        blaze::DynamicVector<double> x0_;
        double tolerance_ = 1e-15;
        double maxIter_ = 10000000;

    public:
        ISparseLinearSolver(const std::vector<std::vector<double>> &A, const std::vector<double> &b);

        virtual ~ISparseLinearSolver() = default;

        virtual std::vector<double> solve() = 0;
    };

    class GaussSeidel : public ISparseLinearSolver {

    public:
        // using base class constructor:
        using ISparseLinearSolver::ISparseLinearSolver;

        std::vector<double> solve() override;
    };

    class BiCGSTAB : public ISparseLinearSolver {

    public:
        // using base class constructor:
        using ISparseLinearSolver::ISparseLinearSolver;

        std::vector<double> solve() override;
    };


}
