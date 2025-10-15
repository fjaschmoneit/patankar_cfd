#include <gtest/gtest.h>
#include "LinEqsSolvers.h"
#include "test_Structs.h"


TEST_F(NK_matrixBuilder, sparseMatrix1_BiCGSTAB)
{
    using SparseMatrix = blaze::CompressedMatrix<LINALG::scalar>;

    unsigned N = 2000;
    SparseMatrix A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    solve_BiCGSTAB<SparseMatrix>( A, b, x, tolerance, maxIter);

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

TEST_F(NK_matrixBuilder, sparseMatrix2_BiCGSTAB)
{
    unsigned N = 2000;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_2<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    solve_BiCGSTAB( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}


TEST_F(NK_matrixBuilder, denseMatrix1_BiCGSTAB)
{
    unsigned N = 600;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setDenseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    solve_BiCGSTAB( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

TEST_F(NK_matrixBuilder, sparseMatrix1_GaussSeidel)
{
    unsigned N = 200;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    solve_GaussSeidel( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

TEST_F(NK_matrixBuilder, sparseMatrix2_GaussSeidel)
{
    unsigned N = 200;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_2<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    solve_GaussSeidel( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}


TEST_F(NK_matrixBuilder, denseMatrix1_GaussSeidel)
{
    unsigned N = 200;
    blaze::CompressedMatrix<LINALG::scalar> A(N,N, 5*N);
    blaze::DynamicVector<LINALG::scalar> b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setDenseProblem_1<blaze::CompressedMatrix<LINALG::scalar>>(A, b, solution);

    solve_GaussSeidel( A, b, x, tolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

