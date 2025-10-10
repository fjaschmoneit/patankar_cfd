#include <numeric>
#include <gtest/gtest.h>
#include "sparseMatrixSolvers.h"
#include "denseMatrixSolvers.h"

#include "TypeDefs_NumericsKernel.h"
#include <blaze/Math.h>

void timer(const std::string& command) {
    static std::chrono::high_resolution_clock::time_point start_time;

    if (command == "start") {
        start_time = std::chrono::high_resolution_clock::now();
    }
    else if (command == "end") {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Time: " << duration.count() << " ms\n";
    }
}

template<typename MatrixType>
void fillBand(blaze::Band<MatrixType> band, LINALG::scalar value) {
    for(size_t i = 0; i < band.size(); i++)  band[i] = value;
}

template<typename MatrixType>
void fillBand(blaze::Band<MatrixType> band, LINALG::vector values) {
    for(size_t i = 0; i < band.size(); i++)  band[i] = values[i];
}

struct NK_matrixBuilder : public::testing::Test {

    LINALG::scalar tolerance = 1e-15;
    unsigned int maxIter = 100000;

    //    ( 2 1 0 0 0 0 )    ( x_0 )       (1)
    //    ( 0 2 1 0 0 0 )    ( x_1 )       (2)
    //    ( 0 0 2 1 0 0 )    ( x_2 )       (1)
    //    ( 0 0 0 2 1 0 )  * ( x_3 )    =  (2)
    //    ( 0 0 0 0 2 1 )    ( x_4 )       (1)
    //    ( 0 0 0 0 0 2 )    ( x_5 )       (2)
    //
    // solution: x = (0, 1, 0, 1, 0, 1)ˆT
    template<typename MatrixType>
    void setSparseProblem_1(MatrixType& A, blaze::DynamicVector<LINALG::scalar>& b, blaze::DynamicVector<LINALG::scalar>& solution ) {

        fillBand<MatrixType>( blaze::band(A,0), 2.0 );
        fillBand<MatrixType>( blaze::band(A,1), 1.0 );

        std::fill(b.begin(), b.end(), 1.0);
        std::fill(solution.begin(), solution.end(), 0.0);
        for (int i = 0; i < b.size(); i++)
        {
            if (i % 2 != 0) {
                b[i]++;
                solution[i]++;
            }
        }
    }


    //    ( -2 0 0 0 1 0 0 0 0 0 0 0 ...)    ( x_0 )       (3)
    //    ( 0 -2 0 0 0 1 0 0 0 0 0 0 ...)    ( x_1 )       (2)
    //    ( 0 0 -2 0 0 0 1 0 0 0 0 0 ...)    ( x_2 )       (1)
    //    ( 0 0 0 -2 0 0 0 1 0 0 0 0 ...)  * ( x_3 )    =  (0)
    //    ( 1 0 0 0 0 -2 0 0 0 1 0 0 ...)    ( x_4 )       (0)
    //    ( 0 1 0 0 0 0 -2 0 0 0 1 0 ...)    ( x_5 )       (0)
    //              ......
    //    ( ... 0 0 0 1 0 0 0 -2 0 0 0 1)    ( x_N-5 )      (0)
    //    ( ... 0 0 0 0 1 0 0 0 -2 0 0 0)    ( x_N-4 )      (-N-1)
    //    ( ... 0 0 0 0 0 1 0 0 0 -2 0 0)    ( x_N-3 )      (-N-2)
    //    ( ... 0 0 0 0 0 0 1 0 0 0 -2 0)    ( x_N-2 )      (-N-3)
    //    ( ... 0 0 0 0 0 0 0 1 0 0 0 -2)    ( x_N-1 )      (-N-4)
    //    solution: x = (1, 2, 3, 4, ... , N)ˆT
    template<typename MatrixType>
    void setSparseProblem_2(MatrixType& A, blaze::DynamicVector<LINALG::scalar>& b, blaze::DynamicVector<LINALG::scalar>& solution )
    {
        //only for problems bigger then 7
        fillBand<MatrixType>( blaze::band(A,0), -2.0 );
        fillBand<MatrixType>( blaze::band(A,4), 1.0 );
        fillBand<MatrixType>( blaze::band(A,-4), 1.0 );

        auto N = b.size();
        b[0] = 3.0;
        b[1] = 2.0;
        b[2] = 1.0;
        b[N-4] = -static_cast<double>(N)-1.0;
        b[N-3] = -static_cast<double>(N)-2.0;
        b[N-2] = -static_cast<double>(N)-3.0;
        b[N-1] = -static_cast<double>(N)-4.0;

        std::iota(solution.begin(), solution.end(), 1.0);
    };

    //    Q = 0.5*N*(N+1)
    //
    //    ( 2 1 1 1 1 ... )    ( x_0 )       ( 1 + Q )
    //    ( 1 2 1 1 1 ... )    ( x_1 )       ( 2 + Q )
    //    ( 1 1 2 1 1 ... )    ( x_2 )   =   ( 3 + Q )
    //    ( 1 1 1 2 1 ... )    ( x_3 )       ( 4 + Q )
    //           ...              ...           ...
    //    ( ... 1 1 1 1 2 )    ( x_(N-1) )   ( N + Q )
    //    solution: x = (1, 2, 3, 4, ... , N)ˆT
    template<typename MatrixType>
    void setDenseProblem_1(MatrixType& A, blaze::DynamicVector<LINALG::scalar>& b, blaze::DynamicVector<LINALG::scalar>& solution )
    {
        auto N = A.rows();
        for (size_t  i = 0; i<N; ++i) {
            for (size_t j = 0; j<N; ++j) {
                A(i,j) = 1.0;
            }
        }
        fillBand<MatrixType>( blaze::band(A,0), 2.0 );

        auto Q = 0.5*N*(N+1);
        std::iota(b.begin(), b.end(), 1+Q);
        std::iota(solution.begin(), solution.end(), 1.0);
    }




};


// ------------------------------- TESTS ------------------------------

TEST_F(NK_matrixBuilder, sparseMatrix1_sparseBiCGSTAB)
{
    unsigned N = 2000;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    Sparse::solve_BiCGSTAB( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], 1e-8);
    }
}

TEST_F(NK_matrixBuilder, sparseMatrix2_sparseBiCGSTAB)
{
    unsigned N = 2000;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_2<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    Sparse::solve_BiCGSTAB( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], 1e-8);
    }
}

TEST_F(NK_matrixBuilder, denseMatrix1_sparseBiCGSTAB)
{
    unsigned N = 200;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setDenseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    Sparse::solve_BiCGSTAB( A, b, x, tolerance, maxIter  );
    //
    // std::cout << A << std::endl;
    // std::cout << b << std::endl;
    // std::cout << x << std::endl;
    // std::cout << solution << std::endl;

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], 1e-8);
    }
}


// TEST_F(NK_dense_N, denseSmallGaussSeidelMatrix)
// {
//     setProblemSize(60);
//     Dense::GaussSeidel linEqs(A,b);
//     auto x = linEqs.solve();
//     for(int i = 0; i < solution.size(); i++)
//     {
//         EXPECT_NEAR(x[i], solution[i],tolerance);
//     }
// }
//
// TEST_F(NK_dense_N, denseSmallJacobiMatrix)
// {
//     setProblemSize(60);
//     Dense::JacobiIter linEqs(A,b);
//     auto x = linEqs.solve();
//     for(int i = 0; i < solution.size(); i++)
//     {
//         EXPECT_NEAR(x[i], solution[i],tolerance);
//     }
// }
//
// TEST_F(NK_sparse_N6, denseMidSizeJacobiMatrix)
// {
//     Dense::JacobiIter linEqs(A,b);
//     auto x = linEqs.solve();
//     EXPECT_EQ(x, solution);
// }
//
// TEST_F(NK_sparse_N, denseMidSizeJacobiMatrix)
// {
//     setProblemSize(100);
//     Dense::JacobiIter linEqs(A,b);
//     auto x = linEqs.solve();
//     for(int i = 0; i < solution.size(); i++)
//     {
//         EXPECT_NEAR(x[i], solution[i],tolerance);
//     }
// }
//
// TEST_F(NK_sparse_N, dense1MidSizeBiCGSTABMatrix)
// {
//     setProblemSize(100);
//     Dense::BiCGSTAB linEqs(A,b);
//     auto x = linEqs.solve();
//     for(int i = 0; i < solution.size(); i++)
//     {
//         EXPECT_NEAR(x[i], solution[i],tolerance);
//     }
// }

//
// TEST_F(NK_sparse_N, sparseMidSizeBiCGSTABMatrix)
// {
//     setProblemSize(100);
//     Sparse::BiCGSTAB linEqs(A,b,4);
//
//     auto x = linEqs.solve();
//     for(int i = 0; i < solution.size(); i++)
//     {
//         EXPECT_NEAR(x[i], solution[i],tolerance);
//     }
// }
//
// TEST_F(NK_sparse_N, CheckSolverCanExecuteTwoTimes)
// {
//
//     setProblemSize(100);
//     Dense::JacobiIter  linEqs_Jacobi            (A,b);
//     Dense::GaussSeidel linEqs_GaussSeidel       (A,b);
//     Dense::BiCGSTAB    linEqs_BiCGSTAB          (A,b);
//     Sparse::BiCGSTAB   linEqs_BiCGSTAB_sparse   (A,b,4);
//
//     auto x_jacobi = linEqs_Jacobi.solve();
//     auto x_gaussSeidel = linEqs_GaussSeidel.solve();
//     auto x_BiCGSTAB = linEqs_BiCGSTAB.solve();
//     auto x_BiCGSTAB_sparse = linEqs_BiCGSTAB_sparse.solve();
//
//     auto x_jacobi_second = linEqs_Jacobi.solve();
//     auto x_gaussSeidel_second = linEqs_GaussSeidel.solve();
//     auto x_BiCGSTAB_second = linEqs_BiCGSTAB.solve();
//     auto x_BiCGSTAB_sparse_second = linEqs_BiCGSTAB_sparse.solve();
//
//     EXPECT_EQ(x_jacobi, x_jacobi_second);
//     EXPECT_EQ(x_gaussSeidel, x_gaussSeidel_second);
//     EXPECT_EQ(x_BiCGSTAB, x_BiCGSTAB_second);
//     EXPECT_EQ(x_BiCGSTAB_sparse, x_BiCGSTAB_sparse_second);
//
//
// }
//
// // why is this test called denseAll? What should this test demonstrate?
// TEST_F(NK_sparse_N, denseAll)
// {
//     setProblemSize(200);
//     Dense::JacobiIter  linEqs_Jacobi          (A,b);
//     Dense::GaussSeidel linEqs_GaussSeidel     (A,b);
//     Dense::BiCGSTAB    linEqs_BiCGSTAB        (A,b);
//     Sparse::BiCGSTAB   linEqs_BiCGSTAB_Sparse (A,b,4);
//
//     timer("start");
//     auto x_jacobi = linEqs_Jacobi.solve();
//     timer("end");
//     timer("start");
//     auto x_gaussSeidel = linEqs_GaussSeidel.solve();
//     timer("end");
//     timer("start");
//     auto x_BiCGSTB = linEqs_BiCGSTAB.solve();
//     timer("end");
//     timer("start");
//     auto x_BiCGSTAB_Sparse = linEqs_BiCGSTAB_Sparse.solve();
//     timer("end");
// }

//
// TEST_F(NK_sparse_NCoefficient, Using_Ap_Aw_An_As_Ae_as_input)
// {
//     auto N = 3;
//     setProblemSize(N);
//
//     Sparse::BiCGSTAB   linEqs_BiCGSTAB_Sparse (N);
//     linEqs_BiCGSTAB_Sparse.setDirectionalFlux(ap,FVM::CardinalDirection::centre);
//     linEqs_BiCGSTAB_Sparse.setDirectionalFlux(ae,FVM::CardinalDirection::east);
//     linEqs_BiCGSTAB_Sparse.setDirectionalFlux(aw,FVM::CardinalDirection::west);
//     linEqs_BiCGSTAB_Sparse.setDirectionalFlux(as,FVM::CardinalDirection::south);
//     linEqs_BiCGSTAB_Sparse.setDirectionalFlux(an,FVM::CardinalDirection::north);
//     linEqs_BiCGSTAB_Sparse.setDirectionalFlux(b,FVM::CardinalDirection::su);
//
//
//     auto x = linEqs_BiCGSTAB_Sparse.solve();
//     for(int i = 0; i < solution.size(); i++)
//     {
//         EXPECT_NEAR(x[i], solution[i],param::tolerance_);
//     }
//     int theEnd = 0;
//
// }