#include <numeric>
#include <gtest/gtest.h>
#include "sparseMatrixSolvers.h"
#include "denseMatrixSolvers.h"
#include <blaze/Math.h>
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

// to do, move this function to another place
blaze::DynamicMatrix<double> toBlaze(const std::vector<std::vector<double>>& v) {
    if( v.empty() ) return blaze::DynamicMatrix<double>();   // tom matrix

    const size_t m = v.size();       // antal rækker
    const size_t n = v[0].size();    // antal kolonner

    blaze::DynamicMatrix<double> A(m, n);

    for (size_t i = 0; i < m; ++i) {
        // evt. tjek at v[i].size() == n, ellers er input ikke rektangulært
        for (size_t j = 0; j < n; ++j) {
            A(i,j) = v[i][j];
        }
    }

    return A;
}

// to do, move this function to another place
blaze::DynamicVector<double> toBlaze(const std::vector<double>& v) {
    blaze::DynamicVector<double> b( v.size() );           // allokér
    std::copy(v.begin(), v.end(), b.begin());             // kopiér
    return b;
}

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

TEST_F(dense_N, denseSmallGaussSeidel) {
    setProblemSize(60);         //
    Dense::GaussSeidel linEqs(A,b);
    auto x = linEqs.solve();
    EXPECT_EQ(x, solution);
}

TEST_F(sparse_N6, denseMidSizeJacobiMatrix) {
    auto A_blaze = toBlaze(A);
    auto b_blaze = toBlaze(b);

    auto solution_blaze = toBlaze(solution);
    Dense::JacobiIter linEqs(A_blaze,b_blaze);
    timer("start");
    auto x = linEqs.solve1();
    timer("end");
    EXPECT_EQ(x, solution_blaze);
}

TEST_F(sparse_N6, denseMidSizeJacobiSum) {
    Dense::JacobiIter linEqs(A,b);
    timer("start");
    auto x = linEqs.solve();
    timer("end");
    EXPECT_EQ(x, solution);
}

TEST_F(sparse_N, denseMidSizeJacobi)
{
    setProblemSize(20);
    Dense::JacobiIter linEqs(A,b);
    timer("start");
    auto x = linEqs.solve();
    timer("end");
    for(int i = 0; i < x.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],1e-10);
    }
}

TEST_F(sparse_N, denseMidSizeJacobiMatrix)
{
    setProblemSize(100);
    auto A_blaze = toBlaze(A);
    auto b_blaze = toBlaze(b);

    auto solution_blaze = toBlaze(solution);
    Dense::JacobiIter linEqs(A_blaze,b_blaze);
    timer("start");
    auto x = linEqs.solve1();
    timer("end");
    for(int i = 0; i < solution_blaze.size(); i++)
    {
        EXPECT_NEAR(x[i], solution[i],1e-10);
    }
}



TEST_F(sparse_N6, solveDenseJacobiLinearSystem)
{
    Dense::JacobiIter linEqs(A,b);
    EXPECT_EQ(linEqs.solve(), solution);
}

TEST_F(sparse_N6, solveGaussSeidelDenseLinearSystem)
{
    Dense::GaussSeidel linEqs(A,b);
    EXPECT_EQ(linEqs.solve(), solution);
}




TEST_F(sparse_N6, solveSparseLinearSystem){

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
