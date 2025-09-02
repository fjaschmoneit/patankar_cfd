#include <numeric>
#include <gtest/gtest.h>
#include "sparseMatrixSolvers.h"
#include "denseMatrixSolvers.h"

struct linearEquationsFixture : public ::testing::Test
{
    //    ( 2 1 0 0 0 0 )    ( x_0 )       (1)
    //    ( 0 2 1 0 0 0 )    ( x_1 )       (2)
    //    ( 0 0 2 1 0 0 )    ( x_2 )       (1)
    //    ( 0 0 0 2 1 0 )  * ( x_3 )    =  (2)
    //    ( 0 0 0 0 2 1 )    ( x_4 )       (1)
    //    ( 0 0 0 0 0 2 )    ( x_5 )       (2)
    //
    // solution: x = (0, 1, 0, 1, 0, 1)ˆT
    // Test Version 2233e3e3r3r

    int N = 0;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> solution;

    void SetUp() override
    {
        N = 6;
        A.assign(N, std::vector<double>(N, 0.0));
        b.assign(N, 1.0);
        solution.assign(N, 0.0);

        for (int i = 0; i < N; i++)
        {
            A[i][i] = 2.0;
            if (i != N-1)
                A[i][i+1] = 1.0;
            if (i % 2 != 0)
                b[i]++;
        }
        solution.resize(N);
        for (int i = 0; i < N; i++)
        {
            if (i%2 != 0)
            {
                solution[i]++;
            }
        }
    }
};


TEST_F(linearEquationsFixture, solveDenseJacobiLinearSystem)
{

    linearEquationsFixture::SetUp();

    Dense::JacobiIter linEqs(A,b);
    std::vector<double> x = linEqs.solve();

    EXPECT_EQ(x, solution);
}

TEST_F(linearEquationsFixture, solveGuissSeidelDenseLinearSystem)
{

    linearEquationsFixture::SetUp();

    Dense::GaussSeidel linEqs(A,b);
    std::vector<double> x = linEqs.solve();

    EXPECT_EQ(x, solution);
}




TEST_F(linearEquationsFixture, solveSparseLinearSystem){

    int N = 10;

    // diagonal and off-diagonals:
    std::vector<double> d0(N, -2.0);
    std::vector<double> d1(N-1, 1.0);

    std::vector<double> b(N, 0.0);
    b[N-1] = -(N+1);

    Sparse::GaussSeidel linEqs(d1,d0,d1,b);
    std::vector<double> x = linEqs.solve();

    // solutions = [1,2,3,4,5,.....]
    std::vector<double> solution(N);
    std::iota(solution.begin(), solution.end(), 1.00);

    EXPECT_EQ(x, solution);
}

