#include <gtest/gtest.h>
#include "../../NumericsKernel/src/LinEqsSolvers.h"
#include "KERNEL_test_Structs.h"
#include "../modules/Common/Util.h"

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
    KERNEL::vector b, x, solution;
    const int nameW = 12;
    const int colW  = 10;
    int timeLimitMs = 1100;
    std::vector<int> sizes;
    int numberOfBands = 0;
    std::vector<KERNEL::smatrix> vA;
    template<typename MatrixType>

    void setSparseProblem_1(MatrixType& A, KERNEL::vector& b, KERNEL::vector& solution, int N)
    {
        solution.resize( N ); solution= 0.0;
        int nx= static_cast<int>(std::sqrt(N));
        auto rows = A.rows();
        auto ap = blaze::band(A,0);
        //Set all main-diagonal coefficients to -0.2
        for (size_t i = 0; i < ap.size(); ++i) {ap[i] = -0.2;}
        for (int i = 0; i < numberOfBands; i++)
        {
            auto as = blaze::band(A,static_cast<int>(i*nx + 1));
            auto an = blaze::band(A,static_cast<int>(-i*nx - 1));
            auto size_of_bands = ap.size()-(i*nx + 1);

            //Set all second-diagonal coefficients to 0.01
            for (size_t i = 0; i < size_of_bands; ++i){as[i]   = 0.01;}
            for (size_t i = 0; i < size_of_bands; ++i){an[i]   = 0.01;}
        }
        //KERNEL::dmatrix AA = A;
        std::iota(solution.begin(), solution.end(), 0.0);
        b = A*solution;
    }

    auto printRow(const std::string& methodName, const std::vector<long long>& timesMs)const
    {
        std::cout << std::left << std::setw(nameW) << (methodName);
        for (size_t i = 0; i < timesMs.size(); ++i) {
            if (i < timesMs.size()) {
                std::ostringstream cell;
                cell << timesMs[i] << " ms";
                std::cout << std::right << std::setw(colW) << cell.str();
            } else {
                // Print empty cells for methods that are not run for all sizes
                std::cout << std::right << std::setw(colW) << "-";
            }
        }
        std::cout << std::endl;
    };

    void preallocateMemory()
    {
        util::timer timer1;
        std::cout<<"Number Of Bands: "<<numberOfBands*2+1<<std::endl;
        timer1.start();
        // Reserve memory for the vector containing all matrices
        vA.reserve(sizes.size());
        for (int N : sizes)
        {
            vA.emplace_back(N, N);
            vA.back().reserve((numberOfBands*2+1) * N);
            for (std::size_t r = 0; r < N; ++r) vA.back().reserve(r, numberOfBands*2+1);
        }
        auto allocateTime = timer1.stop();
        std::cout<<"finish with allocation, allocation time :"<< allocateTime/1000 <<" Sec"<<std::endl;
    }

    template<typename SolverFunc>
    void runTimedSolver(const std::string& label, SolverFunc solver)
    {
        auto timer = util::timer();
        std::vector<long long> timesSolve;
        for (size_t i = 0; i < sizes.size(); ++i)
        {
            const int N = sizes[i];
            // reset all vectors to 0.
            b.resize( N ); b = 0.0;
            x.resize( N ); x = 0.0;
            solution.resize( N );solution = 0.0;
            setSparseProblem_1<KERNEL::smatrix>(vA[i], b, solution, N);
            timer.start();
            // Solve linear system Ax = b
            solver(vA[i], x, b);
            const auto solveTime = timer.stop();
            timesSolve.push_back(solveTime);
            for (int i = 0; i < x.size(); ++i)
            {
                EXPECT_NEAR(x[i], solution[i], 1e-8);
            }
            if (solveTime > timeLimitMs) {
                break;
            }
        }
        printRow(label, timesSolve);
    }
};

TEST_F(NK_matrixBuilder, performance_BiCGSTAB_sparseMatrix1)
{

    auto timer1 = util::timer();
    auto timer2 = util::timer();
    auto timer3 = util::timer();

    const std::vector<int> sizes = {125,250,500,1000,2000,4000,8000,16000,32000};
    for (int i = 0; i<sizes.size(); i++)
    {
        auto N = sizes[i];
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
    auto timer3 = util::timer();


    //sizes = {9,125,250,500,1000,2000,4000,8000,16000,32000,64000, 128000,240000, 480000};
    sizes = {15,125,250,500,1000,2000, 4000,8000,16000,32000,64000,120000};

    for (int i = 2; i<4; i++)
    {
        std::vector<KERNEL::smatrix>().swap( vA );
        for (auto& A : vA) {
            A.clear();   // delete all nonzeros elements
        }
        numberOfBands = i;
        preallocateMemory();
        //check collection time of the largest matrix
        {
            timer3.start();
            auto N = sizes[sizes.size()-1];
            b.resize( N ); b = 0.0; x.resize( N ); x = 0.0; solution.resize( N );solution = 0.0;
            setSparseProblem_1<KERNEL::smatrix>(vA.back(), b, solution, N);
            timer3.stop();
            std::cout << "Collection time of a "<<sizes.back()<<"x"<<sizes.back()<<" matrix: "<<timer3.stop()<<" ms"<<std::endl;
        }
        // Printing Header
        {
            std::cout << std::left << std::setw(nameW) << "Methods";
            for (size_t i = 0; i < sizes.size(); ++i)
            {
                std::cout << std::right << std::setw(colW) << sizes[i];
            }
            std::cout << std::endl;
        }
        //Record execution time of the BiCGSTAB solver
        runTimedSolver("BiCGSTAB",[&](const KERNEL::smatrix& A, KERNEL::vector& x, const KERNEL::vector& b)
        {
            solve_BiCGSTAB<KERNEL::smatrix>(A, x, b, 1e-13, 1000000);
        });

        // Record execution time of the Blaze solver
        runTimedSolver("Blaze solver",[&](const KERNEL::dmatrix& A, KERNEL::vector& x, const KERNEL::vector& b)
        {
            blaze::solve(A, x, b);
        });

        // Record execution time of the Jacobi solver
        runTimedSolver("Jacobi",[&](const KERNEL::smatrix& A, KERNEL::vector& x,const KERNEL::vector& b)
        {
            solve_Jacobi(A, x, b, 1e-13, 1000000);
        });

        //Record execution time of the GaussSeidel solver
        runTimedSolver("GaussSeidel",[&](const KERNEL::smatrix& A, KERNEL::vector& x, KERNEL::vector& b)
        {
            solve_GaussSeidel(A, x, b, 1e-13, 1000000);
        });

    }
}