#include <gtest/gtest.h>
#include "KERNEL.h"
// #include "FVMCalculateCoefficient.h"


struct testtest : public ::testing::Test {
};

TEST_F(testtest, interfaceTest) {

    ObjectRegistry objReg;
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

//
//
// // constructing a square domain with
// struct FVM_laplaceTests : public ::testing::Test {
//     static constexpr double tolerance = 1e-6;
//     unsigned int nx, ny, nbCells;
//     double cellSpacing, faceArea;
//
//
//     // these are to be fetched from a future MESH class
//     std::vector<unsigned int> cellIndices_North;
//     std::vector<unsigned int> cellIndices_South;
//     std::vector<unsigned int> cellIndices_East;
//     std::vector<unsigned int> cellIndices_West;
//
//     std::vector<double> field;
//
//     Sparse::matrix setUp(unsigned int nbX)
//     {
//         if (nbX%2 == 0) {
//             std::cerr << "ERROR: nx must be uneven" << std::endl;
//         }
//         nx = nbX;
//         ny = nbX;
//         nbCells = nx*ny;
//
//         Sparse::matrix coeffMat(nx, ny);
//
//         cellSpacing = 1./static_cast<double>(nx);
//         faceArea = cellSpacing * cellSpacing;
//         field.assign(nbCells,0.0);
//
//         cellIndices_North.resize(nx);
//         cellIndices_South.resize(nx);
//         cellIndices_East.resize(ny);
//         cellIndices_West.resize(ny);
//
//         for (unsigned int i = 0; i < nx; i++) {
//             cellIndices_North[i] = i;
//             cellIndices_South[i] = i+nx*(ny-1);
//             cellIndices_East[i] = i*nx+nx-1;
//             cellIndices_West[i] = i*nx;
//         }
//         return coeffMat;
//     }
//
//     std::vector<double> getHorizontalCenterLineValues( ) const {
//         std::vector<double> midvalues(nx,0.0);
//         for (unsigned int i=0; i < nx; i++) {
//             midvalues[i] = field[ nx/2 + i*nx ];
//         }
//         return midvalues;
//     }
//
//     std::vector<double> getVerticalCenterLineValues(  ) const {
//         std::vector<double> midvalues(ny,0.0);
//         for (unsigned int i=0; i < ny; i++) {
//             midvalues[i] = field[ (ny/2) * nx + i ];
//         }
//         return midvalues;
//     }
//
//     std::vector<double> getCellCenterCoordinates( unsigned int index ) const {
//         double x = ( 0.5 + (double)(index % nx) ) * cellSpacing;
//         double y = ( 0.5 + (double)(index / nx) ) * cellSpacing;
//         return std::vector{x, y};
//     }
// };
//



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
// TEST_F(FVM_laplaceTests, FVM_localDerichletBCs) {
//
//     auto linEqs = setUp(5);

    // auto xpos = [&](size_t index)
    // {
    //     return -dx/2 + static_cast<GLOBAL::scalar>(index+1)*dx;
    // };

        // ---------  solve problem and write solution in 'field' container ------------

    //
    // // we should not save the boundary values in an extra vector.
    // GLOBAL::vector b, vBoundary_Value_North, vBoundary_Value_East, vBoundary_Value_South, vBoundary_Value_West;
    // for (unsigned int i=0; i < nx; i++)
    // {
    //     vBoundary_Value_South.push_back(0);
    //     vBoundary_Value_North.push_back(xpos(i));
    //     vBoundary_Value_West.push_back(0);
    //     vBoundary_Value_East.push_back(xpos(nx-i-1));
    // }

    //
    // std::vector<GLOBAL::scalar> ae,aw,an,as,ap,sp,su;
    // ae.assign(nx*ny,0.0);
    // aw = an = as = as = ap = sp = su = ae;
    //
    // double de,dw,ds,dn,dp;
    // dp = de = dw = ds = dn = 0.0;
    // auto A = 1;
    // for (unsigned int i=0; i < nx*ny; i++)
    // {
    //         dp = de = dw = ds = dn = A / dx;
    //         ae[i] = aw[i] = an[i] = as[i] = de;
    //         sp[i] = su[i] = 0;
    // }

    //
    // auto DerichletNorth = [&](unsigned i) { return getCellCenterCoordinates(i)[0]; };
    // auto DerichletSouth = [&](unsigned i) { return 0.0; };
    // auto DerichletEast = [&](unsigned i) { return getCellCenterCoordinates(i)[1]; };
    // auto DerichletWest = [&](unsigned i) { return 0.0; };
    //
    // std::vector<double> b(nbCells, 0.);
    // std::vector<double> ai(nbCells, 0.0);
    // std::vector<double> ap(nbCells, 0.0);
    //
    // // EAST
    // std::ranges::fill(ai, faceArea/cellSpacing);
    // std::ranges::transform(ap, ai, ai.begin(), std::minus<>{});
    // for (unsigned int i = 0; i < cellIndices_East.size(); i++) {
    //     ai[i] = 0.0;
    //     b[i] -= 2*faceArea/cellSpacing * DerichletEast(i);
    //     ap[i] -= faceArea/cellSpacing;
    // }
    // linEqs.setDirectionalFlux(ai, 1);   // dir=1   --> east

    // // NORTH
    // std::ranges::fill(ai, faceArea/cellSpacing);
    // std::ranges::transform(ap, ai, ai.begin(), std::minus<>{});
    // for (unsigned int i = 0; i < cellIndices_North.size(); i++) {
    //     ai[i] = 0.0;
    //     b[i] -= 2*faceArea/cellSpacing * DerichletNorth(i);
    //     ap[i] -= faceArea/cellSpacing;
    // }
    // linEqs.setDirectionalFlux(ai, 2);   // dir=2   --> north
    //
    // // WEST
    // std::ranges::fill(ai, faceArea/cellSpacing);
    // std::ranges::transform(ap, ai, ai.begin(), std::minus<>{});
    // for (unsigned int i = 0; i < cellIndices_West.size(); i++) {
    //     ai[i] = 0.0;
    //     b[i] -= 2*faceArea/cellSpacing * DerichletWest(i);
    //     ap[i] -= faceArea/cellSpacing;
    // }
    // linEqs.setDirectionalFlux(ai, 3);   // dir=3   --> west
    //
    // // SOUTH
    // std::ranges::fill(ai, faceArea/cellSpacing);
    // std::ranges::transform(ap, ai, ai.begin(), std::minus<>{});
    // for (unsigned int i = 0; i < cellIndices_South.size(); i++) {
    //     ai[i] = 0.0;
    //     b[i] -= 2*faceArea/cellSpacing * DerichletSouth(i);
    //     ap[i] -= faceArea/cellSpacing;
    // }
    // linEqs.setDirectionalFlux(ai, 4);   // dir=4   --> south
    //
    //
    // linEqs.setDirectionalFlux(ai, 0);   // dir=0   --> centre
    //
    // field = linEqs.solve();

    //
    // apply_Boundary(vBoundary_Value_North,getFaceNumberNorth(nx),an,su,sp,A,dx);
    // apply_Boundary(vBoundary_Value_West, getFaceNumberWest(nx) ,aw,su,sp,A,dx);
    // apply_Boundary(vBoundary_Value_East, getFaceNumberEast(nx) ,ae,su,sp,A,dx);
    // apply_Boundary(vBoundary_Value_South,getFaceNumberSouth(nx),as,su,sp,A,dx);
    // for (unsigned int i=0; i < nx*ny; i++)
    // {
    //    ap[i] = -(ae[i] + aw[i] + as[i] + an[i] - sp[i]);
    // }
    //
    // linEqs.setDirectionalFlux(ap,FVM::CardinalDirection::centre);
    // linEqs.setDirectionalFlux(ae,FVM::CardinalDirection::east);
    // linEqs.setDirectionalFlux(aw,FVM::CardinalDirection::west);
    // linEqs.setDirectionalFlux(as,FVM::CardinalDirection::south);
    // linEqs.setDirectionalFlux(an,FVM::CardinalDirection::north);
    // linEqs.setDirectionalFlux(su,FVM::CardinalDirection::su);
    // field = linEqs.solve();
    //
    // auto midHorizontalEquated = getHorizontalCenterLineValues();
    // auto midVerticalEquated = getVerticalCenterLineValues();
    //
    // // theoretical solution:
    // std::vector<double> solution( nx, 0.0 );
    // for (unsigned int i=0; i < nx; i++) {
    //     auto x = 0.5;
    //     auto y = 1.0 - ( 0.5 + i )*dx;
    //     solution[i] = x*y;
    // }
    //
    // for(int i = 0; i < solution.size(); i++)
    // {
    //     EXPECT_NEAR(midVerticalEquated[i], solution[i],tolerance);
    // }
    // for(int i = 0; i < solution.size(); i++)
    // {
    //     EXPECT_NEAR(midHorizontalEquated[solution.size()-i-1], solution[i],tolerance);
//     // }
//     EXPECT_EQ(1,0);
// }




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