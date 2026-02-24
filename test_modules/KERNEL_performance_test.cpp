#include <gtest/gtest.h>
#include "../../NumericsKernel/src/LinEqsSolvers.h"
#include "KERNEL_test_Structs.h"
#include "../modules/Common/Util.h"
#include "KERNEL.h"

using namespace LINEQSOLVERS;

struct NK_DynamikMatrixBuilder : public::testing::Test
{

    //    ( 2 1 0 0 0 0 )    ( x_0 )       (1)
    //    ( 0 2 1 0 0 0 )    ( x_1 )       (2)
    //    ( 0 0 2 1 0 0 )    ( x_2 )       (1)
    //    ( 0 0 0 2 1 0 )  * ( x_3 )    =  (2)
    //    ( 0 0 0 0 2 1 )    ( x_4 )       (1)
    //    ( 0 0 0 0 0 2 )    ( x_5 )       (2)
    //
    // solution: x = (0, 1, 0, 1, 0, 1)ˆT

    std::vector<std::map<std::string,std::vector<long long>>> matrixTime;
    int timeLimitMs = 900;
    std::vector<int> sizes;
    int numberOfBands = 0;
    int numberOffDiagPairs = 0;
    template<typename MatrixType>
    void setSparseProblem_1(MatrixType& A, KERNEL::vector& b, KERNEL::vector& solution, int nx)
    {
        auto rows = A.rows();
        auto ap = blaze::band(A,0);
        //Set all main-diagonal coefficients to -0.2
        for (size_t i = 0; i < ap.size(); ++i) {ap[i] = -0.2;}
        for (int i = 0; i < numberOffDiagPairs; i++)
        {
            auto as = blaze::band(A,static_cast<int>(i*nx + 1));
            auto an = blaze::band(A,static_cast<int>(-i*nx - 1));
            auto size_of_bands = ap.size()-(i*nx + 1);

            //Set all second-diagonal coefficients to 0.01
            for (size_t ii = 0; ii < size_of_bands; ++ii){as[ii]   = 0.01;}
            for (size_t ii = 0; ii < size_of_bands; ++ii){an[ii]   = 0.01;}
        }
        std::iota(solution.begin(), solution.end(), 0.0);
        b = A*solution;
    }

    static void printTimingTable(
        const std::vector<int>& sizes,
        const std::vector<std::map<std::string,
        std::vector<long long>>>& matrixTime,
        int bandIndex,
        int outputWidth = 8,
        std::ostream& os = std::cout)
    {
        // Header
        os << std::left << std::setw(outputWidth+16) << "N of elements";
        for (size_t i = 0; i < sizes.size(); ++i) {
            const long long N = 1LL * sizes[i] * sizes[i];
            os << std::right << std::setw(outputWidth) << N;
        }
        os << "\n";

        for (std::map<std::string, std::vector<long long>>::const_iterator it = matrixTime[bandIndex].begin();
             it != matrixTime[bandIndex].end(); ++it)
        {
            const std::string& method = it->first;
            const std::vector<long long>& times = it->second;

            os << std::left << std::setw(outputWidth+16) << method;

            for (size_t i = 0; i < sizes.size(); ++i) {
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

    template<typename SolverFunc>
    void runTimedSolver(const std::string& label, KERNEL::smatrix& A, KERNEL::vector& b, KERNEL::vector& x, const int nx, int bandIndex, SolverFunc solver)
    {
        if (matrixTime[bandIndex][label+"_solve"].empty() ||
            matrixTime[bandIndex][label+"_solve"].back() < timeLimitMs)
        {
            auto timer = util::timer();
            // reset all vectors to 0.
            b = 0.0; x = 0.0;
            KERNEL::vector solution(nx*nx, 0.0);
            timer.start();
            setSparseProblem_1<KERNEL::smatrix>(A, b, solution, nx);
            auto collectionTime = timer.stop();
            if (label == "BiCGSTAB")
                matrixTime[bandIndex][label+"_fill_A_matrix"].push_back(collectionTime);
            timer.start();
            // Solve linear system Ax = b
            solver(A, x, b);
            auto solverTimer = timer.stop();
            matrixTime[bandIndex][label+"_solve"].push_back(solverTimer);
        }
    }
};

std::vector<int> generateBands(int nx, int bandnumber)
{
    std::vector<int> bands;

    if (bandnumber == 5)
    {
        bands = { -nx, -1, 0, 1, nx };
    }
    else if (bandnumber == 9)
    {
        bands = { -2*nx, -nx, -2, -1, 0, 1, 2, nx, 2*nx };
    }
    return bands;
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
        solve_BiCGSTAB<KERNEL::smatrix>( A, x, b, 1e-13, 10000);
        auto T3 = timer3.stop();
        auto T2 = timer2.stop();
        auto T1 = timer1.stop();

        std::cout << N << ", \t" << T1 << " ms, \t" << T2 << " ms, \t" << T3 << " ms " << std::endl;
        EXPECT_NEAR(x[2], solution[2], 1e-8);
    }
}

TEST_F(NK_DynamikMatrixBuilder, Execution_time_comparison_of_linear_solvers_for_increasing_system_size)
{
    auto timer            = util::timer();
    sizes = {4,11,15,22,31,44,63,90,126,179,252,369,490, 692, 1000, 1400,2000};

    matrixTime.resize(3);
    std::vector<int> numberOfBandsVector;
    for (int bandIndex = 1; bandIndex<3; bandIndex++)
    {
        numberOfBands = bandIndex*4+1;
        numberOfBandsVector.push_back(numberOfBands);
        numberOffDiagPairs = (numberOfBands - 1) / 2;
        std::cout << "Number Of Bands: " << numberOfBands << std::endl;
        for (int ni = 0; ni<sizes.size(); ni++)
        {
            int nx = sizes[ni];
            const std::size_t N = static_cast<std::size_t>(nx) * static_cast<std::size_t>(nx);

            timer.start();
            auto vec =generateBands(nx, numberOfBands);
            auto A = KERNEL::createPreallocatedSparseMatrix(N,nx,generateBands(nx, numberOfBands));
            KERNEL::vector x(N, 0.0);
            KERNEL::vector b(N, 0.0);
            auto allocatedTimer = timer.stop();
            matrixTime[bandIndex]["allocation"].push_back(allocatedTimer);

            //Record execution time of the BiCGSTAB solver
            runTimedSolver("BiCGSTAB",A,b,x,nx,bandIndex,[&](const KERNEL::smatrix& A, KERNEL::vector& x, const KERNEL::vector& b)
            {
                KERNEL::solve(A, x, b, 1e-13, 1000000, KERNEL::BiCGSTAB);
            });

            // Record execution time of the Blaze solver
            runTimedSolver("Blaze_solver",A,b,x,nx,bandIndex,[&](const KERNEL::dmatrix& A, KERNEL::vector& x, const KERNEL::vector& b)
            {
                KERNEL::solve(A, x, b,1e-13,1000000,KERNEL::Blaze_automatic);
            });

            // Record execution time of the Jacobi solver
            runTimedSolver("Jacobi",A,b,x,nx,bandIndex,[&](const KERNEL::dmatrix& A, KERNEL::vector& x,const KERNEL::vector& b)
            {
                KERNEL::solve(A, x, b, 1e-13, 100000, KERNEL::Jacobi);
            });

            //Record execution time of the GaussSeidel solver
            runTimedSolver("GaussSeidel",A,b,x,nx,bandIndex,[&](const KERNEL::dmatrix& A, KERNEL::vector& x, KERNEL::vector& b)
            {
                KERNEL::solve(A, x, b, 1e-13, 1000000, KERNEL::GaussSeidel);
            });

            //printTimingTable(sizes, matrixTime,bandIndex);
        }
        std::cout << "----------------------------------------------------------------------Number of Band "+ std::to_string(numberOfBandsVector[0]) +"---------------------------------------------------------------\n";
        printTimingTable(sizes, matrixTime,bandIndex);
    }
    //std::cout << "--------------------------------------------------------------------------Results-------------------------------------------------------------------\n";
    //std::cout << "----------------------------------------------------------------------Number of Band "+ std::to_string(numberOfBandsVector[0]) +"---------------------------------------------------------------\n";
    //printTimingTable(sizes, matrixTime,1);
    //std::cout << "----------------------------------------------------------------------Number of Band "+ std::to_string(numberOfBandsVector[1]) +"----------------------------------------------------------------\n";
    //printTimingTable(sizes, matrixTime,2);
}