#include <numeric>
#include <gtest/gtest.h>
#include "sparseMatrixSolvers.h"
#include "denseMatrixSolvers.h"
#include <MathUtils.h>

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


struct param
{
    static constexpr double tolerance_ = 1e-7;
};

void fillBandInDenseMatrix(std::vector<std::vector<double>>& A, const std::vector<double>& q, int bandOffset)
{
    int n = A.size();
    for (int i = 0; i < n; ++i)
    {
        int j = i + bandOffset;
        if (j >= 0 && j < n)
        {
            A[i][j] = q[i];
        }
    }
}

//    ( 2 1 0 0 0 0 )    ( x_0 )       (1)
//    ( 0 2 1 0 0 0 )    ( x_1 )       (2)
//    ( 0 0 2 1 0 0 )    ( x_2 )       (1)
//    ( 0 0 0 2 1 0 )  * ( x_3 )    =  (2)
//    ( 0 0 0 0 2 1 )    ( x_4 )       (1)
//    ( 0 0 0 0 0 2 )    ( x_5 )       (2)
//
// solution: x = (0, 1, 0, 1, 0, 1)ˆT
struct sparse_N6 : public ::testing::Test
{
    const int N = 6;

    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> solution;

    sparse_N6()
    : A(N, std::vector<double>(N)),
    b(N, 1.0),
    solution(N,0.0)
    {
        fillBandInDenseMatrix(A, std::vector<double>(N, 2), 0);
        fillBandInDenseMatrix(A, std::vector<double>(N-1, 1), 1);

        for (int i = 0; i < N; i++)
        {
            if (i % 2 != 0) {
                b[i]++;
                solution[i]++;
            }
        }
    }
};


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


// solution: x = (1, 2, 3, 4, ... , N)ˆT
struct sparse_N : public ::testing::Test
{
    //only for problems bigger then 7
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> solution;

    void setProblemSize( const int N) {
        A.assign(N, std::vector<double>(N,0.0));
        b.assign(N,0.0);
        solution.assign(N,0.0);

        fillBandInDenseMatrix(A, std::vector<double>(N, -2.0), 0);
        fillBandInDenseMatrix(A, std::vector<double>(N, 1.0), 4);
        fillBandInDenseMatrix(A, std::vector<double>(N, 1.0), -4);

        b[0] = 3.0;
        b[1] = 2.0;
        b[2] = 1.0;
        b[N-1] = -N-4.0;
        b[N-2] = -N-3.0;
        b[N-3] = -N-2.0;
        b[N-4] = -N-1.0;

        std::iota(solution.begin(), solution.end(), 1.0);
    }
};



//    Q = 0.5*N*(N+1)
//
//    ( 2 1 1 1 1 ... )    ( x_0 )       ( 1 + Q )
//    ( 1 2 1 1 1 ... )    ( x_1 )       ( 2 + Q )
//    ( 1 1 2 1 1 ... )    ( x_2 )   =   ( 3 + Q )
//    ( 1 1 1 2 1 ... )    ( x_3 )       ( 4 + Q )
//           ...              ...           ...
//    ( ... 1 1 1 1 2 )    ( x_(N-1) )   ( N + Q )
//
// solution: x = (1, 2, 3, 4, ... , N)ˆT
struct dense_N : public ::testing::Test {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> solution;

    void setProblemSize(const int N) {

        A.assign(N, std::vector<double>(N,1.0));
        fillBandInDenseMatrix(A, std::vector<double>(N, 2.0), 0);

        b.assign(N,0.0);
        std::iota(b.begin(), b.end(), 1+0.5*N*(N+1) );

        solution.assign(N,0.0);
        std::iota(solution.begin(), solution.end(), 1.0);
    }
};

TEST_F(dense_N, denseSmallGaussSeidelMatrix)
{
    setProblemSize(60);
    Dense::GaussSeidel linEqs(A,b);
    auto x = linEqs.solve();
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],param::tolerance_);
    }
}

TEST_F(dense_N, denseSmallJacobiMatrix)
{
    setProblemSize(60);
    Dense::JacobiIter linEqs(A,b);
    auto x = linEqs.solve();
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],param::tolerance_);
    }
}

TEST_F(sparse_N6, denseMidSizeJacobiMatrix)
{

    Dense::JacobiIter linEqs(A,b);
    auto x = linEqs.solve();
    EXPECT_EQ(x, solution);
}

TEST_F(sparse_N, denseMidSizeJacobiMatrix)
{
    setProblemSize(100);
    Dense::JacobiIter linEqs(A,b);
    auto x = linEqs.solve();
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],param::tolerance_);
    }
}

TEST_F(sparse_N, dense1MidSizeBiCGSTABMatrix)
{
    setProblemSize(100);
    Dense::BiCGSTAB linEqs(A,b);
    auto x = linEqs.solve();
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],param::tolerance_);
    }
}

TEST_F(sparse_N, sparseMidSizeBiCGSTABMatrix)
{
    setProblemSize(100);
    Sparse::BiCGSTAB linEqs(A,b);
    auto x = linEqs.solve();
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],param::tolerance_);
    }
}

TEST_F(sparse_N, CheckSolverCanExecuteTwoTimes)
{

    setProblemSize(100);
    Dense::JacobiIter  linEqs_Jacobi            (A,b);
    Dense::GaussSeidel linEqs_GaussSeidel       (A,b);
    Dense::BiCGSTAB    linEqs_BiCGSTAB          (A,b);
    Sparse::BiCGSTAB   linEqs_BiCGSTAB_sparse   (A,b);

    auto x_jacobi = linEqs_Jacobi.solve();
    auto x_gaussSeidel = linEqs_GaussSeidel.solve();
    auto x_BiCGSTAB = linEqs_BiCGSTAB.solve();
    auto x_BiCGSTAB_sparse = linEqs_BiCGSTAB_sparse.solve();

    auto x_jacobi_second = linEqs_Jacobi.solve();
    auto x_gaussSeidel_second = linEqs_GaussSeidel.solve();
    auto x_BiCGSTAB_second = linEqs_BiCGSTAB.solve();
    auto x_BiCGSTAB_sparse_second = linEqs_BiCGSTAB_sparse.solve();

    EXPECT_EQ(x_jacobi, x_jacobi_second);
    EXPECT_EQ(x_gaussSeidel, x_gaussSeidel_second);
    EXPECT_EQ(x_BiCGSTAB, x_BiCGSTAB_second);
    EXPECT_EQ(x_BiCGSTAB_sparse, x_BiCGSTAB_sparse_second);


}

// why is this test called denseAll? What should this test demonstrate?
TEST_F(sparse_N, denseAll)
{
    setProblemSize(200);
    Dense::JacobiIter  linEqs_Jacobi          (A,b);
    Dense::GaussSeidel linEqs_GaussSeidel     (A,b);
    Dense::BiCGSTAB    linEqs_BiCGSTAB        (A,b);
    Sparse::BiCGSTAB   linEqs_BiCGSTAB_Sparse (A,b);

    timer("start");
    auto x_jacobi = linEqs_Jacobi.solve();
    timer("end");
    timer("start");
    auto x_gaussSeidel = linEqs_GaussSeidel.solve();
    timer("end");
    timer("start");
    auto x_BiCGSTB = linEqs_BiCGSTAB.solve();
    timer("end");
    timer("start");
    auto x_BiCGSTAB_Sparse = linEqs_BiCGSTAB_Sparse.solve();
    timer("end");
}