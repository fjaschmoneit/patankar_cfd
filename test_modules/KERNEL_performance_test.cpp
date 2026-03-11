#include <gtest/gtest.h>
#include <ranges>
#include <fstream>

#include "../modules/Common/Util.h"
#include "KERNEL.h"
#include "KERNEL_test_Structs.h"

struct NK_DynamikMatrixBuilder : public::testing::Test {

    template<typename MatrixType>
    void setSparseProblem_1(MatrixType& A, KERNEL::vector& b, KERNEL::vector& solution)
    {
        std::vector<int> bandIDs = KERNEL::getBandIDs(A);

        for (auto id : bandIDs) {
            KERNEL::fillBand(blaze::band(A,id), 0.1*abs(id) );
        }
        KERNEL::fillBand(blaze::band(A,0), -static_cast<int>(bandIDs.size()) );

        std::iota(solution.begin(), solution.end(), 0.0);
        b = A*solution;
    }

    // FJA return also number of iterations
    // make output easyly readable csv format
    // inlcude also float/double in test
    template<typename MatrixType, typename SolverFunc, typename... SolverArgs>
    long long timeSolver(
        const MatrixType& A,
        KERNEL::vector& x,
        const KERNEL::vector& b,
        const KERNEL::vector& s,
        SolverFunc solver,
        SolverArgs&&... solverArgs){

        auto timer            = util::timer();
        auto testSampleAddress = x.size()%7;

        int nTests = 4;
        long long timeDiff = 0;
        for (auto i : std::views::iota(0,nTests) ) {
            std::ranges::fill(x,0.0);

            timer.start();
            solver(A, x, b, std::forward<SolverArgs>(solverArgs)...);
            timeDiff += timer.stop();
        }

        EXPECT_NEAR(s[testSampleAddress], x[testSampleAddress], 1e-8);

        return timeDiff/nTests;
    }

    static void printTimingTable(
        const std::vector<size_t>& Ns,
        const std::map<const char*,std::vector<long long>>& timeMap,
        std::ostream& os)
    {
        int outputWidth = 8;

        os << std::left << std::setw(outputWidth+16) << "Nb of elements";
        for ( auto N:Ns){
            os << std::right << std::setw(outputWidth) << N;
        }
        os << "\n";

        for (const auto& [method, times] : timeMap){

            os << std::left << std::setw(outputWidth+16) << method;

            for (size_t i = 0; i < Ns.size(); ++i) {
                if (i < times.size()) {
                    std::ostringstream cell;
                    cell << times[i] << " ms";
                    os << std::right << std::setw(outputWidth) << cell.str();
                } else {
                    os << std::right << std::setw(outputWidth) << "-";
                }
            }
            os << "\n";
        }
    }
};

TEST_F(NK_DynamikMatrixBuilder, Execution_time_comparison_of_linear_solvers_for_increasing_system_size)
{
    std::string matrixType = "sparse";

    const std::vector<size_t> Ns = (matrixType=="dense")
        ? std::vector<size_t>{500,1'000, 5'000, 10'000}
        : std::vector<size_t>{100'000, 500'000, 1'000'000};

    std::vector<std::vector<int>> arrayOfBandIDs{
        {0,1,-1},
        {-10,-1,0,1,10},
        {-20,-10,-1,0,1,10,20},
        // {-40,-20,-10,-1,0,1,10,20}
    };

    std::ofstream os( std::string(TEST_SOURCE_DIR) + "/Performance_Results/performanceTestResults_" +matrixType + "_" + util::timer::today()  + ".txt");
    if (!os.is_open())
        throw std::runtime_error("Could not open output file");

    for (auto bandIDs : arrayOfBandIDs) {
        std::map< const char*,std::vector<long long> > timeMap;

        for (auto N : Ns)
        {
            auto timer              = util::timer();
            long long dt            = 0.0;

            // creating the matrix and vectors
            timer.start();
                KERNEL::smatrix A = KERNEL::newTempBandedSMatrix(N, bandIDs, -1.);
                KERNEL::vector x(N, 0.0);
                KERNEL::vector b(N, 0.0);
                KERNEL::vector s(N, 0.0);
            timeMap["allocation sparse"].push_back(timer.stop());

            // filling matrix with values
            timer.start();
                setSparseProblem_1<KERNEL::smatrix>(A, b, s);
            timeMap["fill_A_matrix sparse"].push_back(timer.stop());

            if (matrixType=="sparse") {
                // I need this function pointer to tell the compiler which solve() overload to choose.
                using SparseSolverPtr = void(*)(const KERNEL::smatrix&,KERNEL::vector&,const KERNEL::vector&,KERNEL::SolverMethod,GLOBAL::scalar,unsigned int);

                dt = timeSolver(A, x, b, s, static_cast<SparseSolverPtr>(KERNEL::solve), KERNEL::BiCGSTAB, 1e-13, 1'000'000u);
                timeMap["BiCGSTAB"].push_back(dt);

            }else {
                const KERNEL::dmatrix Ad(A);
                using DenseSolverPtr = void(*)(const KERNEL::dmatrix&,KERNEL::vector&,const KERNEL::vector&,KERNEL::SolverMethod,GLOBAL::scalar,unsigned int);

                dt = timeSolver(Ad, x, b, s, static_cast<DenseSolverPtr>(KERNEL::solve), KERNEL::BiCGSTAB, 1e-13, 1'000'000u);
                timeMap["BiCGSTAB"].push_back(dt);

                dt = timeSolver(Ad, x, b, s, static_cast<DenseSolverPtr>(KERNEL::solve), KERNEL::Jacobi, 1e-13, 1'000'000u);
                timeMap["Jacobi"].push_back(dt);

                // sth. wrong with gauss seidel. It is very slow.
                // dt = timeSolver(Ad, x, b, s, static_cast<DenseSolverPtr>(KERNEL::solve), KERNEL::GaussSeidel, 1e-13, 1'000'000u);
                // timeMap["Gauss-Seidel"].push_back(dt);

                dt = timeSolver(Ad, x, b, s, static_cast<DenseSolverPtr>(KERNEL::solve), KERNEL::Blaze_automatic, 1e-13, 1'000'000u);
                timeMap["Blaze native"].push_back(dt);

            }
    }

    os << "----------------------------------------------------------------------Bands {";
    std::ranges::for_each( bandIDs, [&](const auto& n) { os << std::scientific << std::setprecision(3)<< ' ' << n; });
    os << " }---------------------------------------------------------------\n";
    printTimingTable(Ns, timeMap, os);
    }
}


TEST_F(NK_matrixBuilder, performance_BiCGSTAB_sparseMatrix1)
{
    auto timer1 = util::timer();
    auto timer2 = util::timer();
    auto timer3 = util::timer();

    const std::vector<int> sizes = {11,15,22,31,44,63,90,126,178};
    for (int i = 0; i<sizes.size(); i++)
    {
        auto nx = sizes[i];
        auto N = nx*nx;
        timer1.start();
        KERNEL::smatrix A(N,N, 5*N);
        KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

        timer2.start();
        setSparseProblem_1<KERNEL::smatrix>(A, b, solution);
        //setSparseProblem_2<KERNEL::smatrix>(A, b, solution);
        //setDenseProblem_1<KERNEL::smatrix>(A, b, solution);

        timer3.start();
        KERNEL::solve(A, x, b, KERNEL::BiCGSTAB);
        // solve_BiCGSTAB<KERNEL::smatrix>( A, x, b, 1e-13, 10000);
        auto T3 = timer3.stop();
        auto T2 = timer2.stop();
        auto T1 = timer1.stop();

        std::cout << N << ", \t" << T1 << " ms, \t" << T2 << " ms, \t" << T3 << " ms " << std::endl;
        EXPECT_NEAR(x[2], solution[2], 1e-8);
    }
}