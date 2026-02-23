//
// Created by Peter Berg Ammundsen on 20/01/2026.
//
#include <gtest/gtest.h>
#include "KERNEL.h"
#include <numeric>
#include "Util.h"
#include "structured2d.h"
#include <iostream>



struct MESHInterface : public ::testing::Test
{
    GLOBAL::scalar lengthX  = 0;
    GLOBAL::scalar lengthY  = 0;
    GLOBAL::scalar firstXCenterPosition  = 0;
    GLOBAL::scalar secondXCenterPosition = 0;
    GLOBAL::scalar firstYCenterPosition  = 0;
    GLOBAL::scalar secondYCenterPosition = 0;
    GLOBAL::scalar dx = 0;
    GLOBAL::scalar dy = 0;

    GLOBAL::scalar LastXCenterPosition = 0;
    GLOBAL::scalar LastYCenterPosition = 0;

    GLOBAL::scalar cellThickness = 0.0;
    GLOBAL::scalar faceAreaX = 0.0;
    GLOBAL::scalar faceAreaY = 0.0;
    void setMeshProblem(unsigned int nbcellsX , unsigned int nbcellsY )
    {
        dx = lengthX/nbcellsX;
        dy = lengthY/nbcellsY;

        firstXCenterPosition  = (0.5*dx);
        secondXCenterPosition = (0.5*dx)+dx;
        firstYCenterPosition  = (nbcellsY-0.5)*dy;
        secondYCenterPosition = (nbcellsY-0.5)*dy;

        LastXCenterPosition = lengthX-0.5*dx;
        LastYCenterPosition = 0.5*dy;

        cellThickness = 1.0;
        faceAreaX = dx*cellThickness;
        faceAreaY = dy*cellThickness;
    }

};

// Validates all public function in the structured2dRegularRectangle mesh class
// for a 2D structured rectangular FVM mesh.
// Layout of the cell near (0,0) with coordinate.
//        y ^
//          |
//         1dy ----- n ------ NE
//          |        |        |
//          |        |        |
//          0.5dy -- P ------ e      → x (East)
//          |        |        |
//          |        |        |
//        --0 ----- 0.5dx---- dx--> x
TEST_F(MESHInterface, interfaceTest1)
{
    lengthX = 1712.0;
    lengthY = 123.0;

    unsigned int nbcellsX = 3312;
    unsigned int nbcellsY = 312;

    setMeshProblem(nbcellsX, nbcellsY);

    MESH::structured2dRegularRectangle mesh(lengthX,nbcellsX,lengthY,nbcellsY);
    mesh.init();
    EXPECT_EQ(mesh.lenX(),lengthX);
    EXPECT_EQ(mesh.lenY(),lengthY);

    EXPECT_EQ(mesh.nbCellsX(),nbcellsX);
    EXPECT_EQ(mesh.nbCellsY(),nbcellsY);
    EXPECT_EQ(mesh.nbCells() ,nbcellsX*nbcellsY);

    EXPECT_EQ(mesh.getCellCenterCoordinate_X(0) , firstXCenterPosition);
    EXPECT_EQ(mesh.getCellCenterCoordinate_X(1) , secondXCenterPosition);
    EXPECT_EQ(mesh.getCellCenterCoordinate_Y(0) , firstYCenterPosition);
    EXPECT_EQ(mesh.getCellCenterCoordinate_Y(1) , secondYCenterPosition);

    EXPECT_NEAR(mesh.getCellCenterCoordinate_X(mesh.nbCells()-1) , LastXCenterPosition,1e-10);
    EXPECT_NEAR(mesh.getCellCenterCoordinate_Y(mesh.nbCells()-1) , LastYCenterPosition,1e-10);

    EXPECT_EQ(mesh.getCellSpacing_X() , dx);
    EXPECT_EQ(mesh.getCellSpacing_Y() , dy);

    EXPECT_EQ(mesh.getCellFaceArea_X() , faceAreaX);
    EXPECT_EQ(mesh.getCellFaceArea_Y() , faceAreaY);

    EXPECT_NEAR(mesh.getCellReciprocalSpacing_X() , 1/dx , 1e-10);
    EXPECT_NEAR(mesh.getCellReciprocalSpacing_Y() , 1/dy   , 1e-10);

    //checking boundary cellID.
    std::vector<int >   NorthBoundaryCellNumber(nbcellsX),SouthBoundaryCellNumber(nbcellsX),
                         WestBoundaryCellNumber(nbcellsX), EastBoundaryCellNumber(nbcellsX);
    std::iota(NorthBoundaryCellNumber.begin(),NorthBoundaryCellNumber.end(),0);
    WestBoundaryCellNumber = SouthBoundaryCellNumber = EastBoundaryCellNumber = NorthBoundaryCellNumber;
    std::ranges::transform(SouthBoundaryCellNumber,SouthBoundaryCellNumber.begin(),[nbcellsX,nbcellsY] (auto v) { return v+(nbcellsY-1)*nbcellsX ; });
    std::ranges::transform(WestBoundaryCellNumber ,WestBoundaryCellNumber.begin() ,[nbcellsX]          (auto v) { return v*nbcellsX;               });
    std::ranges::transform(EastBoundaryCellNumber ,EastBoundaryCellNumber.begin() ,[nbcellsX]          (auto v) { return (v+1)*nbcellsX -1 ;       });

    //checking boundary top cells
    {
        for (auto cellId : mesh.region(MESH::RegionID::Boundary_top))
        {
            bool isCellIdFound = std::ranges::find(NorthBoundaryCellNumber, cellId) != NorthBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
        }
        //checking getCellFacePos function
        for (const auto& [cellId, pos] : mesh.getCellFacesPos(MESH::RegionID::Boundary_top))
        {
            GLOBAL::scalar x = cellId * dx + 0.5*dx;
            GLOBAL::scalar y = lengthY;
            bool isCellIdFound = std::ranges::find(NorthBoundaryCellNumber, cellId) != NorthBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
            EXPECT_DOUBLE_EQ(pos.x,x);
            EXPECT_DOUBLE_EQ(pos.y,y);
        }
    }
    //checking boundary bottom cells
    {
        for (auto cellId : mesh.region(MESH::RegionID::Boundary_bottom))
        {
            bool isCellIdFound = std::ranges::find(SouthBoundaryCellNumber, cellId) != SouthBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
        }
        //checking getCellFacePos function
        for (const auto& [cellId, pos] : mesh.getCellFacesPos(MESH::RegionID::Boundary_bottom))
        {
            GLOBAL::scalar x =   (cellId - (nbcellsX*nbcellsY) + nbcellsX) * dx + 0.5*dx;
            GLOBAL::scalar y = 0;
            bool isCellIdFound = std::ranges::find(SouthBoundaryCellNumber, cellId) != SouthBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
            EXPECT_DOUBLE_EQ(pos.x,x);
            EXPECT_DOUBLE_EQ(pos.y,y);
        }
    }
    //checking boundary left cells
    {
        for (auto cellId : mesh.region(MESH::RegionID::Boundary_left))
        {
            bool isCellIdFound = std::ranges::find(WestBoundaryCellNumber, cellId) != WestBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
        }
        //checking getCellFacePos function
        for (const auto& [cellId, pos] : mesh.getCellFacesPos(MESH::RegionID::Boundary_left))
        {
            GLOBAL::scalar x = 0;
            GLOBAL::scalar y = lengthY-0.5*dy-dy*cellId / nbcellsX;
            bool isCellIdFound = std::ranges::find(WestBoundaryCellNumber, cellId) != WestBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
            EXPECT_FLOAT_EQ(pos.x , x);
            EXPECT_FLOAT_EQ(pos.y , y);
        }
    }
    //checking boundary right cells
    {
        for (auto cellId : mesh.region(MESH::RegionID::Boundary_right))
        {
            bool isCellIdFound = std::ranges::find(EastBoundaryCellNumber, cellId) != EastBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
        }
        //checking getCellFacePos function
        for (const auto& [cellId, pos] : mesh.getCellFacesPos(MESH::RegionID::Boundary_right))
        {
            GLOBAL::scalar x =  lengthX;
            GLOBAL::scalar y = lengthY-0.5*dy-dy*((cellId-nbcellsX+1)/nbcellsX);
            bool isCellIdFound = std::ranges::find(EastBoundaryCellNumber, cellId) != EastBoundaryCellNumber.end();
            EXPECT_TRUE(isCellIdFound);
            EXPECT_TRUE(mesh.isBoundaryCell(cellId));
            EXPECT_FLOAT_EQ(pos.x , x);
            EXPECT_FLOAT_EQ(pos.y ,y);
        }
    }
}