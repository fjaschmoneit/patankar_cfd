#include <gtest/gtest.h>
#include "KERNEL.h"
#include <numeric>

#include "../../NumericsKernel/src/LinEqsSolvers.h"
#include "blaze/Blaze.h"
// #include "FVMCalculateCoefficient.h"


struct kernelInterface : public ::testing::Test {
};

TEST_F(kernelInterface, interfaceTest) {

    KERNEL::ObjectRegistry objReg;
    auto vecHandle = objReg.newVector(5, 3);
    auto matHandleDense = objReg.newMatrix(5, 3, false);
    auto matHandleSparse = objReg.newMatrix(5, 3, true);

    // I have to close the registry, before I am allowed to access its objects:
    EXPECT_THROW(objReg.getVectorRef(vecHandle),std::runtime_error);

    // now it is closed.
    objReg.closeRegistry();

    // These are references. They may go out of scope, the objects survive in the reg.
    auto u = objReg.getVectorRef(vecHandle);
    auto A  = objReg.getDenseMatrixRef(matHandleDense);
    auto B  = objReg.getSparseMatrixRef(matHandleSparse);

    // It is forbidden to add new objects in the registry, for keeping the existing
    // references valid. The following line will fail:
    EXPECT_THROW(objReg.newMatrix(5, 3, false),std::runtime_error);

}

//
//
// std::vector<unsigned int> getFaceNumberSouth(size_t size_of_vector)
// {
//     std::vector<unsigned int> values(size_of_vector);
//     std::iota(values.begin(), values.end(), static_cast<int>(size_of_vector * (size_of_vector-1)));
//     return values;
// }
// std::vector<unsigned int> getFaceNumberNorth(size_t size_of_vector)
// {
//     std::vector<unsigned int> values(size_of_vector);
//     std::iota(values.begin(), values.end(), static_cast<int>(0));
//     return values;
// }
// std::vector<unsigned int> getFaceNumberWest(size_t size_of_vector)
// {
//     std::vector<unsigned int> values(size_of_vector);
//     std::iota(values.begin(), values.end(), static_cast<int>(0));
//     for (auto &x : values | std::views::all)
//     {
//         x *= size_of_vector;
//     }
//     return values;
// }
//
// std::vector<unsigned int> getFaceNumberEast(size_t size_of_vector)
// {
//     std::vector<unsigned int> values(size_of_vector);
//     std::iota(values.begin(), values.end(), static_cast<int>(0));
//     for (auto &x : values | std::views::all)
//     {
//         x = x * (size_of_vector)+ size_of_vector-1;
//     }
//     return values;
// }
// face ændre navn til bounadry cells indicies
//
// void apply_Boundary(const std::vector<FVM::scalar>& value, std::vector<unsigned int> Faces, std::vector<FVM::scalar>& input, std::vector<FVM::scalar>& Su, std::vector<FVM::scalar>& Sp,double A,double dx)
// {
//     unsigned int i = 0;
//     //apply Direchtly boundary condition.
//     for (const unsigned int face : Faces)
//     {
//         Sp[face] += -2 * A / dx;
//         Su[face] +=  -2 * A / dx * value[i];
//         input[face] = 0;
//         i++;
//     }
// }



// constructing a square domain with
struct FVM_laplaceTests : public ::testing::Test {
    static constexpr double tolerance = 1e-6;
    unsigned int nx, ny, nbCells;
    double lenx, leny, cellSpacing, faceArea;
    KERNEL::MatrixHandle AHandle;
    KERNEL::VectorHandle uHandle;
    KERNEL::VectorHandle bHandle;

    // these are to be fetched from a future MESH class
    std::vector<unsigned int> cellIndices_North;
    std::vector<unsigned int> cellIndices_South;
    std::vector<unsigned int> cellIndices_East;
    std::vector<unsigned int> cellIndices_West;

    KERNEL::ObjectRegistry setUp(KERNEL::scalar length, unsigned int nbX)
    {
        if (nbX%2 == 0)
            throw std::runtime_error("ERROR: nx must be uneven");

        nx = nbX;
        ny = nbX;
        nbCells = nx*ny;
        leny = length;
        lenx = length;
        KERNEL::ObjectRegistry objReg;
        AHandle = objReg.newMatrix(nbCells, nbCells, true);
        uHandle = objReg.newVector(nbCells);
        bHandle = objReg.newVector(nbCells);
        objReg.closeRegistry();

        cellSpacing = lenx/static_cast<double>(nx);
        faceArea = cellSpacing * cellSpacing;

        cellIndices_North.resize(nx);
        cellIndices_South.resize(nx);
        cellIndices_East.resize(ny);
        cellIndices_West.resize(ny);

        for (unsigned int i = 0; i < nx; i++) {
            cellIndices_North[i] = i;
            cellIndices_South[i] = i+nx*(ny-1);
            cellIndices_East[i] = i*nx+nx-1;
            cellIndices_West[i] = i*nx;
        }
        return objReg;
    }

    // std::vector<double> getHorizontalCenterLineValues( ) const {
    //     std::vector<double> midvalues(nx,0.0);
    //     for (unsigned int i=0; i < nx; i++) {
    //         midvalues[i] = field[ nx/2 + i*nx ];
    //     }
    //     return midvalues;
    // }
    //
    std::vector<double> getVerticalCenterLineValues( KERNEL::vector& field ) const {
        std::vector<double> midvalues(ny,0.0);
        for (unsigned int i=0; i < ny; i++) {
            midvalues[i] = field[ (ny/2) * nx + i ];
        }
        return midvalues;
    }

};


TEST_F(FVM_laplaceTests, FVM_testtest) {

    // nb of cells along one side
    auto nbX = 81;

    // setup() fills the registry with a matrix and a vector,
    // the handles are members of the test struct
    auto objReg = setUp(1.0, nbX);

    auto u = objReg.getVectorRef(uHandle);
    auto b = objReg.getVectorRef(bHandle);
    auto A = objReg.getDenseMatrixRef(AHandle);
    // auto A = objReg.getSparseMatrixRef(AHandle);

    // building matrix
    auto diagonal = blaze::band(A,0);
    for (unsigned int i = 0; i < diagonal.size(); i++) {
        diagonal[i] = 4.0;
    }

    auto aw = blaze::band(A,-1);
    for (unsigned int i = 0; i < aw.size(); i++) {
        aw[i] = 1.0;
    }

    auto ae = blaze::band(A,1);
    for (unsigned int i = 0; i < ae.size(); i++) {
        ae[i] = 2.0;
    }

    auto correctResult = *KERNEL::newTempVector(b.size());
    std::iota(correctResult.begin(),correctResult.end(),1);

    b = A*correctResult;

    // an educated guess for u:
    std::ranges::fill(u, 2.0);

    KERNEL::solve(A, u, b, 1e-10, 1000, KERNEL::BiCGSTAB);

    for (unsigned int i = 0; i < correctResult.size(); i++) {
        EXPECT_NEAR(u[i], correctResult[i], 1e-4);
    }
}



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
TEST_F(FVM_laplaceTests, FVM_localDerichletBCs) {

    auto objReg = setUp(1.0, 11);

    auto A = objReg.getSparseMatrixRef(AHandle);
    auto u = objReg.getVectorRef(uHandle);
    auto b = objReg.getVectorRef(bHandle);

    KERNEL::vector ae(nbCells, 0.0),aw(nbCells, 0.0),an(nbCells, 0.0),as(nbCells, 0.0),ap(nbCells, 0.0),sp(nbCells, 0.0),su(nbCells, 0.0);

    for (unsigned int i=0; i < nbCells; i++)
    {
            ae[i] = aw[i] = an[i] = as[i] = faceArea / cellSpacing;
            ap[i] = -4 * faceArea / cellSpacing;
    }

    // EAST
    for (unsigned int i = 0; i < ny; i++) {
        auto xpos = lenx;
        auto ypos = leny - (0.5+i)*cellSpacing;
        auto j = nx-1 +i*nx;
        ae[j] = 0.0;
        b[j] -= 2*faceArea/cellSpacing * ypos*xpos;
        ap[j] -= faceArea/cellSpacing;
    }

    // NORTH
    for (unsigned int i = 0; i < nx; i++) {
        auto xpos = (0.5*cellSpacing + i*cellSpacing);
        auto ypos = leny;
        an[i] = 0.0;
        b[i] -= 2*faceArea/cellSpacing * xpos*ypos;
        ap[i] -= faceArea/cellSpacing;
    }

    // WEST
    for (unsigned int i = 0; i < ny; i++) {
        auto xpos = 0.0;
        auto ypos = (leny - 0.5*cellSpacing - i * cellSpacing);
        auto j = i*nx;
        aw[j] = 0.0;
        b[j] -= 2*faceArea/cellSpacing * xpos*ypos;
        ap[j] -= faceArea/cellSpacing;
    }

    // SOUTH
    for (unsigned int i = 0; i < nx; i++) {
        auto xpos = (0.5*cellSpacing + i*cellSpacing);
        auto ypos = 0.0;
        auto j = (ny-1)*nx + i;
        as[j] = 0.0;
        b[j] -= 2*faceArea/cellSpacing * xpos*ypos;
        ap[j] -= faceArea/cellSpacing;
    }

    for (unsigned int i = 0; i < ap.size(); i++) {
        blaze::band(A,0)[i] = ap[i];
    }
    for (unsigned int i = 0; i < ap.size() -1; i++) {
        blaze::band(A,1)[i] = ae[i];
    }
    for (unsigned int i = 0; i < ap.size()-1; i++) {
        blaze::band(A,-1)[i] = aw[i+1];
    }
    for (unsigned int i = 0; i < ap.size()-nx; i++) {
        blaze::band(A,nx)[i] = as[i];
    }
    int mnx = -1*nx;        // I don't see the necessity for that
    for (unsigned int i = 0; i < ap.size()-nx; i++) {
        blaze::band(A,mnx)[i] = an[i+nx];
    }

    KERNEL::solve(A, u, b, 1e-15, 1000, KERNEL::BiCGSTAB);

    // theoretical solution, vertical mid-line at x = lenx/2
    KERNEL::vector solution( nx, 0.0 );
    for (unsigned int i=0; i < nx; i++) {
        auto x = 0.5*lenx;
        auto y = leny - ( 0.5 + i )*cellSpacing;
        solution[i] = x*y;
    }

    for(int i = 0; i < solution.size(); i++)
    {
        auto j = nx/2 + nx*i;
        EXPECT_NEAR(u[j], solution[i],tolerance);
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
//TEST_F(FVM_laplaceTests, spacVarDerichletBCs) {
    // setUp(51);
    //
    // // ---------  solve problem and write solution in 'field' container ------------
    // auto xpos = [&](size_t index)
    // {
    //     return -dx/2 + static_cast<GLOBAL::scalar>(index+1)*dx;
    // };
    // GLOBAL::vector b, vBoundary_Value_North, vBoundary_Value_East, vBoundary_Value_South, vBoundary_Value_West;
    // vBoundary_Value_South.resize(nx);
    // vBoundary_Value_North = vBoundary_Value_West = vBoundary_Value_East = vBoundary_Value_South;
    // for (unsigned int i=0; i < nx; i++)
    // {
    //     vBoundary_Value_South[i] = xpos(i)*xpos(i);
    //     vBoundary_Value_North[i] = xpos(i)*xpos(i)-1;
    //     vBoundary_Value_West[nx-i-1]  = -xpos(i)*xpos(i);
    //     vBoundary_Value_East[nx-i-1]  = 1-xpos(i)*xpos(i);
    // }
    //
    // GLOBAL::vector ae,aw,an,as,ap,sp,su;
    // ae.assign(nx*ny,0.0);
    // aw = an = as = as = ap = sp = su = ae;
    //
    // double de,dw,ds,dn,dp;
    // dp = de = dw = ds = dn = 0.0;
    // auto A = 1;
    // for (unsigned int i=0; i < nx*ny; i++)
    // {
    //     dp = de = dw = ds = dn = A / dx;
    //     ae[i] = aw[i] = an[i] = as[i] = de;
    //     sp[i] = su[i] = 0;
    // }
    //
    // Sparse::BiCGSTAB linEqs(static_cast<int>(nx));
    //
    // apply_Boundary(vBoundary_Value_East, getFaceNumberEast(nx) ,ae,su,sp,A,dx);
    // linEqs.setDirectionalFlux(ae,FVM::CardinalDirection::east);
    //
    // apply_Boundary(vBoundary_Value_West, getFaceNumberWest(nx) ,aw,su,sp,A,dx);
    // linEqs.setDirectionalFlux(aw,FVM::CardinalDirection::west);
    //
    // apply_Boundary(vBoundary_Value_South,getFaceNumberSouth(nx),as,su,sp,A,dx);
    // linEqs.setDirectionalFlux(as,FVM::CardinalDirection::south);
    //
    // apply_Boundary(vBoundary_Value_North,getFaceNumberNorth(nx),an,su,sp,A,dx);
    // linEqs.setDirectionalFlux(an,FVM::CardinalDirection::north);
    // linEqs.setDirectionalFlux(su,FVM::CardinalDirection::su);
    //
    // for (unsigned int i=0; i < nx*ny; i++)
    // {
    //     ap[i] = -(ae[i] + aw[i] + as[i] + an[i] - sp[i]);
    // }
    // linEqs.setDirectionalFlux(ap,FVM::CardinalDirection::centre);
    //
    // field = linEqs.solve();
    //
    // auto midHorizontalEquated = getHorizontalCenterLineValues();
    // auto midVerticalEquated = getVerticalCenterLineValues();
    //
    // // theoretical solution:
    // std::vector<double> solution( nx, 0.0 );
    // for (unsigned int i=0; i < nx; i++) {
    //     auto x = 0.5;
    //     auto y = ( 0.5 + i )*dx;
    //     solution[i] = x*x - y*y;
    // }
    //
    // for(int i = 0; i < solution.size(); i++)
    // {
    //     EXPECT_NEAR(midVerticalEquated[i], solution[i],tolerance);
    // }
    // for(int i = 0; i < solution.size(); i++)
    // {
    //     EXPECT_NEAR(midHorizontalEquated[i], solution[solution.size()-i-1],tolerance);
    // }
// }


// 2D Poisson Equation Test Case (Fitzpatrick Example)
//
// Problem Setup:
// - Domain: 0 ≤ x ≤ 1, 0 ≤ y ≤ 1 (rectangular)
// - PDE: ∇²φ = f(x,y)
// - Source term: f(x,y) = 6xy (1-y) - 2x^3
//
// Boundary Conditions:
// - Dirichlet: φ(x,0) = φ(x,1) = 0, φ(0,y) = 0,  φ(1,y) = y(1-y)
//
// Analytical Solution:
// φ(x,y) = y*(1-y)x^3
//
// Reference: R. Fitzpatrick, "An example solution of Poisson's equation in 2-d"
// https://farside.ph.utexas.edu/teaching/329/lectures/node71.html
TEST_F(FVM_laplaceTests, 2DPoissonDerichlet) {


}