#include <gtest/gtest.h>
#include "../../NumericsKernel/src/LinEqsSolvers.h"
#include "KERNEL_test_Structs.h"

using namespace LINEQSOLVERS;

// This test validates the BiCGSTAB solver on a simple sparse upper-bidiagonal system (main diagonal 2, superdiagonal 1)
// where the exact solution is known. The RHS is constructed accordingly, and we check that the computed
// solution matches the analytical one within tolerance
TEST_F(NK_matrixBuilder, sparseMatrix1_BiCGSTAB)
{
    unsigned N = 2000;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<KERNEL::smatrix>(A, b, solution);

    solve_BiCGSTAB<KERNEL::smatrix>( A, x, b, AlgoTolerance, maxIter);

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], TestTolerance);
    }
}

// This test validates the BiCGSTAB solver on a sparse banded system with a 3-band stencil (main diagonal value at -2 and ±4 off-diagonals, value of 1).
// The RHS is constructed so the analytical solution is known (x = [1,2,3,…,N]), and we verify the computed solution matches within tolerance.
TEST_F(NK_matrixBuilder, sparseMatrix2_BiCGSTAB)
{
    unsigned N = 2000;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_2<KERNEL::smatrix>(A, b, solution);

    solve_BiCGSTAB( A, x, b, AlgoTolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], TestTolerance);
}

// This test validates the BiCGSTAB solver on a dense linear system with 2 on the diagonal and 1 elsewhere.
// The right-hand side is constructed so the exact solution is known (x = [1,2,3,…,N]) and verified within tolerance.
TEST_F(NK_matrixBuilder, denseMatrix1_BiCGSTAB)
{
    unsigned N = 600;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setDenseProblem_1<KERNEL::smatrix>(A, b, solution);

    solve_BiCGSTAB( A, x, b, AlgoTolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], TestTolerance);
}

// This test validates the Jacobi solver on a simple sparse upper-bidiagonal system (diag=2, superdiag=1) with a known exact solution.
// We solve Ax=b from a zero initial guess and verify the computed solution matches the analytical one within tolerance.
TEST_F(NK_matrixBuilder, sparseMatrix1_Jacobi)
{
    unsigned N = 200;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<KERNEL::smatrix>(A, b, solution);

    solve_Jacobi( A, x, b, AlgoTolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], TestTolerance);
    }
}

// This test validates the Jacobi solver on the sparse 3-band system with main diagonal (value -2) and ±4 off-diagonals (value 1).
// We solve Ax=b for a problem with known exact solution (x = [1,2,3,…,N]) and verify the computed result within tolerance.
TEST_F(NK_matrixBuilder, sparseMatrix2_Jacobi)
{
    unsigned N = 200;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_2<KERNEL::smatrix>(A, b, solution);

    solve_Jacobi( A, x, b, AlgoTolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], TestTolerance);
}

// This test checks that the Jacobi method does NOT converge for the dense system from setDenseProblem_1 (diag=2, off-diagonals=1).
// We do not expect Jacobi to converge for this dense system, because the iteration matrix has spectral radius ρ > 1.
// Therefore doesJacobiConverge(A) should return false and the solve step is skipped.
TEST_F(NK_matrixBuilder, denseMatrix1_Jacobi)
{
    unsigned N = 40;
    KERNEL::dmatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);
    setDenseProblem_1<KERNEL::dmatrix>(A, b, solution);
    bool doesConverge = doesJacobiConverge(A);
    EXPECT_FALSE(doesConverge);
    if(doesConverge)
    {
        solve_Jacobi( A,x, b, AlgoTolerance, maxIter  );
        for(int i = 0; i < solution.size(); i++)
        {
            EXPECT_NEAR(x[i], solution[i], TestTolerance);
        }
    }else
    {
        std::cout<<"Skipping solving, as the Jacobi method does not converge."<<std::endl;
    }

}

// TEST_F(NK_matrixBuilder, sparseMatrix1_GaussSeidel)
// {
//     unsigned N = 200;
//     KERNEL::smatrix A(N,N, 5*N);
//     KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);
//
//     setSparseProblem_1<KERNEL::smatrix>(A, b, solution);
//
//     solve_GaussSeidel( A, x, b, tolerance, maxIter  );
//
//     for(int i = 0; i < solution.size(); i++)
//         EXPECT_NEAR(x[i], solution[i], 1e-8);
// }

// TEST_F(NK_matrixBuilder, sparseMatrix2_GaussSeidel)
// {
//     unsigned N = 200;
//     KERNEL::smatrix A(N,N, 5*N);
//     KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);
//
//     setSparseProblem_2<KERNEL::smatrix>(A, b, solution);
//
//     solve_GaussSeidel( A, x, b, tolerance, maxIter  );
//
//     for(int i = 0; i < solution.size(); i++)
//         EXPECT_NEAR(x[i], solution[i], 1e-8);
// }

// This test validates the Gauss–Seidel solver on the dense system from setDenseProblem_1 (diag=2, off-diagonals=1).
// We solve Ax=b for a case with known exact solution (x = [1,2,3,…,N]) and verify the result within tolerance.
TEST_F(NK_matrixBuilder, denseMatrix1_GaussSeidel)
{
    unsigned N = 40;
    KERNEL::dmatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setDenseProblem_1<KERNEL::dmatrix>(A, b, solution);

    solve_GaussSeidel( A, x, b, AlgoTolerance, maxIter  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], TestTolerance);
}

