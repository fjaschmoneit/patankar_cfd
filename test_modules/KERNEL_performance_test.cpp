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
    const int outputWidth = 12;
    int timeLimitMs = 900;
    std::vector<int> sizes;
    int numberOfBands = 0;
    int numberOffDiagPairs = 0;

    KERNEL::MatrixHandle AHandle;
    KERNEL::VectorHandle uHandle;
    KERNEL::VectorHandle bHandle;

    KERNEL::ObjectRegistry  AllocatedMatrixAndVector(unsigned int nx)
    {
        const int N = static_cast<int>(nx*nx);
        KERNEL::ObjectRegistry objReg;

        AHandle = objReg.newMatrix(N, N, true);
        uHandle = objReg.newVector(N);
        bHandle = objReg.newVector(N);

        objReg.closeRegistry();

        auto& A = objReg.getSparseMatrixRef(AHandle);
        A.resize(N, N);
        A.reset();

        return objReg;
    }
    template<typename MatrixType>
    void setSparseProblem_1(MatrixType& A, KERNEL::vector& b, KERNEL::vector& solution, int nx)
    {
        int N = nx*nx;
        A.reset();
        solution.resize( N ); solution= 0.0;
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
        int outputWidth = 12,
        std::ostream& os = std::cout)
    {
        // Header
        os << std::left << std::setw(outputWidth) << "N of elements";
        for (size_t i = 0; i < sizes.size(); ++i) {
            const long long N = 1LL * sizes[i] * sizes[i];
            os << std::right << std::setw(outputWidth - 2) << N;
        }
        os << "\n";

        for (std::map<std::string, std::vector<long long>>::const_iterator it = matrixTime[bandIndex].begin();
             it != matrixTime[bandIndex].end(); ++it)
        {
            const std::string& method = it->first;
            const std::vector<long long>& times = it->second;

            os << std::left << std::setw(outputWidth) << method;

            for (size_t i = 0; i < sizes.size(); ++i) {
                if (i < times.size()) {
                    std::ostringstream cell;
                    cell << times[i] << " ms";
                    os << std::right << std::setw(outputWidth - 2) << cell.str();
                } else {
                    os << std::right << std::setw(outputWidth - 2) << "-";
                }
            }
            os << "\n";
        }
    }

    void preallocateMemory(KERNEL::smatrix& A,
                           std::size_t N,
                           std::size_t& reservedRows)
    {
        A.reserve(static_cast<std::size_t>(numberOfBands) * N);

        for (std::size_t r = reservedRows; r < N; ++r)
        {
            A.reserve(r, numberOfBands);
        }
        reservedRows = N;
    }

    template<typename SolverFunc>
    void runTimedSolver(const std::string& label, KERNEL::smatrix& A, KERNEL::vector& b, KERNEL::vector& x, const int nx, int bandIndex, SolverFunc solver)
    {
        if (matrixTime[bandIndex][label].empty() ||
            matrixTime[bandIndex][label].back() < timeLimitMs)
        {
            auto timer = util::timer();
            const int N = nx*nx;
            // reset all vectors to 0.
            b.resize( N ); b = 0.0;
            x.resize( N ); x = 0.0;
            KERNEL::vector solution;
            solution.resize( N );solution = 0.0;
            setSparseProblem_1<KERNEL::smatrix>(A, b, solution, nx);
            timer.start();
            // Solve linear system Ax = b
            solver(A, x, b);
            matrixTime[bandIndex][label].push_back(timer.stop());
        }
    }
};

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
    sizes = {4,11,15,22,31,44,63,90,126,179,252,369,490, 692};
    matrixTime.resize(3);
    for (int bandIndex = 1; bandIndex<3; bandIndex++)
    {
        numberOfBands = bandIndex*4+1;
        numberOffDiagPairs = (numberOfBands - 1) / 2;
        std::cout << "Number Of Bands: " << numberOfBands << std::endl;

        auto objReg = AllocatedMatrixAndVector(static_cast<unsigned int>(sizes.back()));
        auto& x = objReg.getVectorRef(uHandle);
        auto& b = objReg.getVectorRef(bHandle);
        auto& A = objReg.getSparseMatrixRef(AHandle);

        // Track how much has been reserved so far (monotonic growth)
        std::size_t reservedRows  = 0;
        std::size_t reservedTotal = 0;

        for (int ni = 0; ni<sizes.size(); ni++)
        {
            int nx = sizes[ni];
            const std::size_t N = static_cast<std::size_t>(nx) * static_cast<std::size_t>(nx);

            timer.start();
            A.resize(N, N);
            A.reset();
            b.resize(N); b = 0.0;
            x.resize(N); x = 0.0;
            preallocateMemory(A, N, reservedRows);
            matrixTime[bandIndex]["allocation"].push_back(timer.stop());

            //Record execution time of the BiCGSTAB solver
            runTimedSolver("BiCGSTAB",A,b,x,nx,bandIndex,[&](const KERNEL::smatrix& A, KERNEL::vector& x, const KERNEL::vector& b)
            {
                KERNEL::solve(A, x, b, 1e-13, 1000000, KERNEL::BiCGSTAB);
            });

            // Record execution time of the Blaze solver
            runTimedSolver("Blaze solver",A,b,x,nx,bandIndex,[&](const KERNEL::dmatrix& A, KERNEL::vector& x, const KERNEL::vector& b)
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

            printTimingTable(sizes, matrixTime,bandIndex, outputWidth);
        }
        std::cout << "--------------------------------------------------------------------------Results-------------------------------------------------------------------\n";
        std::cout << "----------------------------------------------------------------------Number of Band 5---------------------------------------------------------------\n";
        printTimingTable(sizes, matrixTime,1, outputWidth);
        std::cout << "----------------------------------------------------------------------Number of Band 9----------------------------------------------------------------\n";
        printTimingTable(sizes, matrixTime,2, outputWidth);
    }
}