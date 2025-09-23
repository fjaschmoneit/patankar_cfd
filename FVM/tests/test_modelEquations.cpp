#include <valarray>
#include <gtest/gtest.h>
#include "denseMatrixSolvers.h"

struct laplaceTests : public ::testing::Test {
    unsigned int nx=5;
    unsigned int ny=5;
    double dx=0.2;

    std::vector<double> field;

    void setUp(unsigned int size)
    {
        if (nx%2 == 0) {
            std::cerr << "ERROR: nx must be uneven" << std::endl;
        }
        nx = size;
        ny = size;
        dx = 1./static_cast<double>(nx);
        field.assign(nx*ny,0.0);
    }

    std::vector<double> getHorizontalCenterLineValues( ) const {
        std::vector<double> midvalues(nx,0.0);
        for (unsigned int i=0; i < nx; i++) {
            midvalues[i] = field[ nx/2 + i*nx ];
        }
        return midvalues;
    }

    std::vector<double> getVerticalCenterLineValues(  ) const {
        std::vector<double> midvalues(ny,0.0);
        for (unsigned int i=0; i < ny; i++) {
            midvalues[i] = field[ (ny/2) * nx + i ];
        }
        return midvalues;
    }

    std::vector<double> getCellCenterCoordinates( unsigned int index ) const {
        double x = ( 0.5 + (double)(index % nx) ) * dx;
        double y = ( 0.5 + (double)(index / nx) ) * dx;
        return std::vector{x, y};
    }
};




// 2D Laplace Equation with constant BCs
//
// Problem Setup:
// - Domain: 0 ≤ x ≤ 1, 0 ≤ y ≤ 1 (rectangular)
// - ODE: ∇²φ = 0
//
// Boundary Conditions:
// - Dirichlet: φ(x,0) = 0, φ(0,y) = 0, φ(x,1) = 1,  φ(1,y) = 1
//
// Analytical Solution:
// φ(x,y) = y*x
TEST_F(laplaceTests, constDerichletBCs) {
    setUp(21);

    // ---------  solve problem and write solution in 'field' container ------------




    // -------------------------------------------------------------------

    auto midHorizontalEquated = getHorizontalCenterLineValues();
    auto midVerticalEquated = getVerticalCenterLineValues();

    // theoretical solution:
    std::vector<double> solution( nx, 0.0 );
    for (unsigned int i=0; i < nx; i++) {
        auto x = 0.5;
        auto y = ( 0.5 + i )*dx;
        solution[i] = x*y;
    }

    EXPECT_EQ(midHorizontalEquated, solution);
    EXPECT_EQ(midVerticalEquated, solution);
}



// 2D Laplace Equation with spacially varying BCs
//
// Problem Setup:
// - Domain: 0 ≤ x ≤ 1, 0 ≤ y ≤ 1 (rectangular)
// - ODE: ∇²φ = 0
//
// Boundary Conditions:
// - Dirichlet: φ(x,0) = xˆ2 , φ(0,y) = -yˆ2, φ(x,1) = xˆ2-1,  φ(1,y) = 1-yˆ2
//
// Analytical Solution:
// φ(x,y) = xˆ2-yˆ2
TEST_F(laplaceTests, spacVarDerichletBCs) {
    setUp(21);

    // ---------  solve problem and write solution in 'field' container ------------




    // -------------------------------------------------------------------

    auto midHorizontalEquated = getHorizontalCenterLineValues();
    auto midVerticalEquated = getVerticalCenterLineValues();

    // theoretical solution:
    std::vector<double> solution( nx, 0.0 );
    for (unsigned int i=0; i < nx; i++) {
        auto x = 0.5;
        auto y = ( 0.5 + i )*dx;
        solution[i] = x*x - y*y;
    }

    EXPECT_EQ(midVerticalEquated, solution);
}