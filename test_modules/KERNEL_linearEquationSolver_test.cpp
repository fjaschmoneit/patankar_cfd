#include <gtest/gtest.h>
#include "../../NumericsKernel/src/LinEqsSolvers.h"
#include "KERNEL.h"


#include "KERNEL_test_Structs.h"
using namespace LINEQSOLVERS;

TEST_F(NK_matrixBuilder, sparseMatrix1_BiCGSTAB)
{
    unsigned N = 2000;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<KERNEL::smatrix>(A, b, solution);

    solve_BiCGSTAB<KERNEL::smatrix>( A, x, b, 1e-13, 10000);

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], 1e-8);
    }
}

// This test validates the BiCGSTAB solver on a sparse banded system with a 3-band stencil (main diagonal value at -2 and ±4 off-diagonals, value of 1).
// The RHS is constructed so the analytical solution is known (x = [1,2,3,…,N]), and we verify the computed solution matches within tolerance.
TEST_F(NK_matrixBuilder, sparseMatrix2_BiCGSTAB)
{
    unsigned N = 200;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_2<KERNEL::smatrix>(A, b, solution);

    solve_BiCGSTAB( A, x, b, 1e-15, 10000  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

// This test validates the BiCGSTAB solver on a dense linear system with 2 on the diagonal and 1 elsewhere.
// The right-hand side is constructed so the exact solution is known (x = [1,2,3,…,N]) and verified within tolerance.
TEST_F(NK_matrixBuilder, denseMatrix1_BiCGSTAB)
{
    unsigned N = 600;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setDenseProblem_1<KERNEL::smatrix>(A, b, solution);

    solve_BiCGSTAB( A, x, b, 1e-14, 10000  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

// This test validates the Jacobi solver on a simple sparse upper-bidiagonal system (diag=2, superdiag=1) with a known exact solution.
// We solve Ax=b from a zero initial guess and verify the computed solution matches the analytical one within tolerance.
TEST_F(NK_matrixBuilder, sparseMatrix1_Jacobi)
{
    unsigned N = 200;
    KERNEL::smatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);

    setSparseProblem_1<KERNEL::smatrix>(A, b, solution);

    solve_Jacobi( A, x, b, 1e-15, 10000  );

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i], 1e-8);
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

    solve_Jacobi( A, x, b, 1e-15, 20000  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

// This test checks that the Jacobi method does NOT converge for the dense system from setDenseProblem_1 (diag=2, off-diagonals=1).
// We do not expect Jacobi to converge for this dense system, because the iteration matrix has spectral radius ρ > 1.
// Therefore solve_Jacobi should return false and the check solution step is skipped.
TEST_F(NK_matrixBuilder, denseMatrix1_Jacobi)
{
    unsigned N = 40;
    KERNEL::dmatrix A(N,N, 5*N);
    KERNEL::vector b(N,0.0), x(N, 0.0), solution(N, 0.0);
    setDenseProblem_1<KERNEL::dmatrix>(A, b, solution);
    bool jacobiConverged = true;
    try {
        solve_Jacobi( A,x, b, 1e-15, 10000  );
    } catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        jacobiConverged = false;
    }

    EXPECT_FALSE(jacobiConverged); //
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

    solve_GaussSeidel( A, x, b, 1e-15, 10000  );

    for(int i = 0; i < solution.size(); i++)
        EXPECT_NEAR(x[i], solution[i], 1e-8);
}

