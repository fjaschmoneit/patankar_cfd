#include <gtest/gtest.h>
#include "LinEqsSolvers.h"
#include "test_Structs.h"

class timer {

    private:

    std::chrono::high_resolution_clock::time_point start_time;

    public:

    void start(){
        start_time = std::chrono::high_resolution_clock::now();
    }

    long long stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    }
};


TEST_F(NK_matrixBuilder, performance_BiCGSTAB_sparseMatrix1)
{
    using SparseMatrix = blaze::CompressedMatrix<LINALG::scalar>;
    auto timer1 = timer();
    auto timer2 = timer();
    auto timer3 = timer();

    std::vector<int> sizes = {1000,2000,4000,8000,16000,32000,64000};
    for (int i = 0; i<sizes.size(); i++) {
        auto N = sizes[i];
        timer1.start();
        SparseMatrix A(N,N, 5*N);
        blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

        timer2.start();
        setSparseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);
        //setSparseProblem_2<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);
        //setDenseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

        timer3.start();
        solve_BiCGSTAB<SparseMatrix>( A, b, x, tolerance, maxIter);
        auto T3 = timer3.stop();
        auto T2 = timer2.stop();
        auto T1 = timer1.stop();

        std::cout << N << ", \t" << T1 << " ms, \t" << T2 << " ms, \t" << T3 << " ms " << std::endl;
        EXPECT_NEAR(x[2], solution[2], 1e-8);
    }
}


// // why is this test called denseAll? What should this test demonstrate?
 TEST_F(NK_matrixBuilder, Sparse_run_all_tests)
 {
    auto timer1 = timer();
    auto timer2 = timer();
    auto timer3 = timer();

    using SparseMatrix = blaze::CompressedMatrix<LINALG::scalar>;
    int size = 1000;
    auto N = size;
    SparseMatrix A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);
    setSparseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    timer1.start();
    solve_GaussSeidel<SparseMatrix>( A, b, x, tolerance, maxIter);
    auto T1 = timer1.stop();
    timer2.start();
    solve_JacobiIter<SparseMatrix>( A, b, x, tolerance, maxIter);
    auto T2 = timer2.stop();
    timer3.start();
    solve_BiCGSTAB<SparseMatrix>( A, b, x, tolerance, maxIter);
    auto T3 = timer3.stop();

    std::cout << N << ", \t" << T1 << " ms, \t" << T2 << " ms, \t" << T3 << " ms " << std::endl;
 }

//TEST_F(NK_matrixBuilder, Dense_run_all_tests)
//{
//    auto timer1 = timer();
//    auto timer2 = timer();
//    auto timer3 = timer();
//
//    using DenseMatrix = blaze::DynamicMatrix<LINALG::scalar>;
//    int size = 1000;
//    auto N = size;
//    DenseMatrix A(N,N, 5*N);
//    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);
//    setSparseProblem_1<blaze::DynamicMatrix<LINALG::scalar>>(A, b, solution);
//
//    timer1.start();
//    //solve_GaussSeidel<DenseMatrix>( A, b, x, tolerance, maxIter);
//    auto T1 = timer1.stop();
//    timer2.start();
//    solve_BiCGSTAB( A, b, x, tolerance, maxIter);
//    auto T2 = timer2.stop();
//    timer3.start();
//    solve_BiCGSTAB<DenseMatrix>( A, b, x, tolerance, maxIter);
//    auto T3 = timer3.stop();
//
//    std::cout << N << ", \t" << T1 << " ms, \t" << T2 << " ms, \t" << T3 << " ms " << std::endl;
//
//
//}
