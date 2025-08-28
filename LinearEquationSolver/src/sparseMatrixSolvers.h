
#include <vector>

namespace Sparse{

    // interface
    class ISparseLinearSolver {

    protected:
        std::vector<std::vector<double>> A_;
        std::vector<double> b_;

    public:


        ISparseLinearSolver(std::vector<double> aw, std::vector<double> ap, std::vector<double> ae, std::vector<double> b)
        : b_(b) {
            // creating sparse matrix with a_i diagonals
        }


        virtual ~ISparseLinearSolver() = default;

        virtual std::vector<double> solve() = 0;
    };

    class GaussSeidel : public ISparseLinearSolver {

    public:
        // using base class constructor:
        using ISparseLinearSolver::ISparseLinearSolver;

        std::vector<double> solve() override;
    };


}
