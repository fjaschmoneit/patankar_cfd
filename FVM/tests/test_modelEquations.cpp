#include <valarray>
#include <gtest/gtest.h>
#include "denseMatrixSolvers.h"
#include "sparseMatrixSolvers.h"
#include "FVMCalculateCoefficient.h"

struct param
{
    static constexpr double tolerance_ = 1e-4;
};

std::vector<unsigned int> getFaceNumberSouth(size_t size_of_vector)
{
    std::vector<unsigned int> values(size_of_vector);
    std::iota(values.begin(), values.end(), static_cast<int>(size_of_vector * (size_of_vector-1)));
    return values;
}
std::vector<unsigned int> getFaceNumberNorth(size_t size_of_vector)
{
    std::vector<unsigned int> values(size_of_vector);
    std::iota(values.begin(), values.end(), static_cast<int>(0));
    return values;
}
std::vector<unsigned int> getFaceNumberWest(size_t size_of_vector)
{
    std::vector<unsigned int> values(size_of_vector);
    std::iota(values.begin(), values.end(), static_cast<int>(0));
    for (auto &x : values | std::views::all)
    {
        x *= size_of_vector;
    }
    return values;
}

std::vector<unsigned int> getFaceNumberEast(size_t size_of_vector)
{
    std::vector<unsigned int> values(size_of_vector);
    std::iota(values.begin(), values.end(), static_cast<int>(0));
    for (auto &x : values | std::views::all)
    {
        x = x * (size_of_vector)+ size_of_vector-1;
    }
    return values;
}
// face ændre navn til bounadry cells indicies
void apply_Boundary(const GLOBAL::vector& value,std::vector<unsigned int> Faces,GLOBAL::vector& input,GLOBAL::vector& Su,GLOBAL::vector& Sp,double A,double dx)
{
    unsigned int i = 0;
    //apply Direchtly boundary condition.
    for (const unsigned int face : Faces)
    {
        Sp[face] += -2 * A / dx;
        Su[face] +=  -2 * A / dx * value[i];
        input[face] = 0;
        i++;
    }
}

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




// 2D Laplace Equation with local BCs
//
// Problem Setup:
// - Domain: 0 ≤ x ≤ 1, 0 ≤ y ≤ 1 (rectangular)
// - ODE: ∇²φ = 0
//
// Boundary Conditions:
// - Dirichlet: φ(x,0) = 0, φ(0,y) = 0, φ(x,1) = x,  φ(1,y) = y
//
// Analytical Solution:
// φ(x,y) = y*x
TEST_F(laplaceTests, localDerichletBCs) {
    //setUp(21);
    setUp(41);
    auto L = 1;
    auto xpos = [&](size_t index)
    {
        return -dx/2 + static_cast<GLOBAL::scalar>(index+1)*dx;
    };
    /*
        // ---------  solve problem and write solution in 'field' container ------------


        auto L = 0.5;
        auto xpos = [&](size_t index)
        {
            auto dx = L / nx;
            return -dx/2 + static_cast<GLOBAL::scalar>(index+1)*dx;
        };
        GLOBAL::vector ai;
        GLOBAL::vector b, vBoundary_Value_East_And_North, vBoundary_Value_South_And_West;

        for (unsigned int i=0; i < nx; i++)
        {
            //vBoundary_Value_West_And_North.push_back(xpos(i));
            vBoundary_Value_East_And_North.push_back(1);
            vBoundary_Value_South_And_West.push_back(0);
        }

        FVMCalculateCoefficient fvmCoef(nx,ny,nx*ny,static_cast<GLOBAL::scalar>(L));

        fvmCoef.calculateDiffusionAidx(FVM::CardinalDirection::east);
        fvmCoef.calculateDiffusionAidx(FVM::CardinalDirection::west);
        fvmCoef.calculateDiffusionAidy(FVM::CardinalDirection::south);
        fvmCoef.calculateDiffusionAidy(FVM::CardinalDirection::north);

        fvmCoef.apply_Boundary(FVM::CardinalDirection::east ,vBoundary_Value_East_And_North);
        fvmCoef.apply_Boundary(FVM::CardinalDirection::north,vBoundary_Value_East_And_North);

        fvmCoef.apply_Boundary(FVM::CardinalDirection::west ,vBoundary_Value_South_And_West);
        fvmCoef.apply_Boundary(FVM::CardinalDirection::south,vBoundary_Value_South_And_West);

        fvmCoef.CollectAp();

        auto ap = fvmCoef.getAi(FVM::CardinalDirection::centre);
        auto ae = fvmCoef.getAi(FVM::CardinalDirection::east);
        auto aw = fvmCoef.getAi(FVM::CardinalDirection::west);
        auto as = fvmCoef.getAi(FVM::CardinalDirection::south);
        auto an = fvmCoef.getAi(FVM::CardinalDirection::north);
        b = fvmCoef.getB();


        Sparse::BiCGSTAB linEqs(b,static_cast<int>(nx));

        linEqs.setDirectionalFlux(ap,FVM::CardinalDirection::centre);
        linEqs.setDirectionalFlux(ae,FVM::CardinalDirection::east);
        linEqs.setDirectionalFlux(aw,FVM::CardinalDirection::west);
        linEqs.setDirectionalFlux(as,FVM::CardinalDirection::south);
        linEqs.setDirectionalFlux(an,FVM::CardinalDirection::north);

        field = linEqs.solve();
    */
    GLOBAL::vector b, vBoundary_Value_North, vBoundary_Value_East, vBoundary_Value_South, vBoundary_Value_West;
    for (unsigned int i=0; i < nx; i++)
    {
        vBoundary_Value_South.push_back(0);
        vBoundary_Value_North.push_back(xpos(i));
        vBoundary_Value_West.push_back(0);
        vBoundary_Value_East.push_back(xpos(nx-i-1));
    }

    GLOBAL::vector ae,aw,an,as,ap,sp,su;
    ae.assign(nx*ny,0.0);
    aw = an = as = as = ap = sp = su = ae;

    double de,dw,ds,dn,dp;
    dp = de = dw = ds = dn = 0.0;
    auto A = 1;
    for (unsigned int i=0; i < nx*ny; i++)
    {
            dp = de = dw = ds = dn = A / dx;
            ae[i] = aw[i] = an[i] = as[i] = de;
            sp[i] = su[i] = 0;
    }


    apply_Boundary(vBoundary_Value_North,getFaceNumberNorth(nx),an,su,sp,A,dx);
    apply_Boundary(vBoundary_Value_West, getFaceNumberWest(nx) ,aw,su,sp,A,dx);
    apply_Boundary(vBoundary_Value_East, getFaceNumberEast(nx) ,ae,su,sp,A,dx);
    apply_Boundary(vBoundary_Value_South,getFaceNumberSouth(nx),as,su,sp,A,dx);
    for (unsigned int i=0; i < nx*ny; i++)
    {
       ap[i] = -(ae[i] + aw[i] + as[i] + an[i] - sp[i]);
    }

    Sparse::BiCGSTAB linEqs(static_cast<int>(nx));
    linEqs.setDirectionalFlux(ap,FVM::CardinalDirection::centre);
    linEqs.setDirectionalFlux(ae,FVM::CardinalDirection::east);
    linEqs.setDirectionalFlux(aw,FVM::CardinalDirection::west);
    linEqs.setDirectionalFlux(as,FVM::CardinalDirection::south);
    linEqs.setDirectionalFlux(an,FVM::CardinalDirection::north);
    linEqs.setDirectionalFlux(su,FVM::CardinalDirection::su);
    field = linEqs.solve();

    auto midHorizontalEquated = getHorizontalCenterLineValues();
    auto midVerticalEquated = getVerticalCenterLineValues();

    // theoretical solution:
    std::vector<double> solution( nx, 0.0 );
    for (unsigned int i=0; i < nx; i++) {
        auto x = 0.5;
        auto y = 1.0 - ( 0.5 + i )*dx;
        solution[i] = x*y;
    }

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(midVerticalEquated[i], solution[i],param::tolerance_);
    }
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(midHorizontalEquated[solution.size()-i-1], solution[i],param::tolerance_);
    }
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
    //setUp(21);
    setUp(51);
    //setUp(3);

    // ---------  solve problem and write solution in 'field' container ------------
    auto xpos = [&](size_t index)
    {
        return -dx/2 + static_cast<GLOBAL::scalar>(index+1)*dx;
    };
    GLOBAL::vector b, vBoundary_Value_North, vBoundary_Value_East, vBoundary_Value_South, vBoundary_Value_West;
    vBoundary_Value_South.resize(nx);
    vBoundary_Value_North = vBoundary_Value_West = vBoundary_Value_East = vBoundary_Value_South;
    for (unsigned int i=0; i < nx; i++)
    {
        vBoundary_Value_South[i] = xpos(i)*xpos(i);
        vBoundary_Value_North[i] = xpos(i)*xpos(i)-1;
        vBoundary_Value_West[nx-i-1]  = -xpos(i)*xpos(i);
        vBoundary_Value_East[nx-i-1]  = 1-xpos(i)*xpos(i);
    }

    GLOBAL::vector ae,aw,an,as,ap,sp,su;
    ae.assign(nx*ny,0.0);
    aw = an = as = as = ap = sp = su = ae;

    double de,dw,ds,dn,dp;
    dp = de = dw = ds = dn = 0.0;
    auto A = 1;
    for (unsigned int i=0; i < nx*ny; i++)
    {
        dp = de = dw = ds = dn = A / dx;
        ae[i] = aw[i] = an[i] = as[i] = de;
        sp[i] = su[i] = 0;
    }

    Sparse::BiCGSTAB linEqs(static_cast<int>(nx));

    apply_Boundary(vBoundary_Value_East, getFaceNumberEast(nx) ,ae,su,sp,A,dx);
    linEqs.setDirectionalFlux(ae,FVM::CardinalDirection::east);

    apply_Boundary(vBoundary_Value_West, getFaceNumberWest(nx) ,aw,su,sp,A,dx);
    linEqs.setDirectionalFlux(aw,FVM::CardinalDirection::west);

    apply_Boundary(vBoundary_Value_South,getFaceNumberSouth(nx),as,su,sp,A,dx);
    linEqs.setDirectionalFlux(as,FVM::CardinalDirection::south);

    apply_Boundary(vBoundary_Value_North,getFaceNumberNorth(nx),an,su,sp,A,dx);
    linEqs.setDirectionalFlux(an,FVM::CardinalDirection::north);
    linEqs.setDirectionalFlux(su,FVM::CardinalDirection::su);

    for (unsigned int i=0; i < nx*ny; i++)
    {
        ap[i] = -(ae[i] + aw[i] + as[i] + an[i] - sp[i]);
    }
    linEqs.setDirectionalFlux(ap,FVM::CardinalDirection::centre);

    field = linEqs.solve();

    auto midHorizontalEquated = getHorizontalCenterLineValues();
    auto midVerticalEquated = getVerticalCenterLineValues();

    // theoretical solution:
    std::vector<double> solution( nx, 0.0 );
    for (unsigned int i=0; i < nx; i++) {
        auto x = 0.5;
        auto y = ( 0.5 + i )*dx;
        solution[i] = x*x - y*y;
    }

    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(midVerticalEquated[i], -solution[i],param::tolerance_);
    }
    for(int i = 0; i < solution.size(); i++)
    {
        EXPECT_NEAR(midHorizontalEquated[i], solution[solution.size()-i-1],param::tolerance_);
    }
}