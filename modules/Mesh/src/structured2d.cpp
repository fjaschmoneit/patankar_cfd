//
// Created by Peter Berg Ammundsen on 20/01/2026.
//
#include "../include/structured2d.h"

using namespace MESH;

// constructor
structured2dRegularRectangle::structured2dRegularRectangle(GLOBAL::scalar lengthX, unsigned int nbCellsX,GLOBAL::scalar lengthY, unsigned int nbCellsY)
: BaseMesh(lengthX, nbCellsX, lengthY, nbCellsY)
{
    if (nbCellsX < 2 || nbCellsY < 2)
    {
        throw std::runtime_error("Mesh dimensions must be >= 2");
    }
    init();
}
void structured2dRegularRectangle::init()
{
    // Initialize boundaries container with 4 sides for a rectangle: 0=bottom,1=left,2=top,3=right
    fillRegion( RegionID::Entire );
    fillRegion( RegionID::Boundary_bottom );
    fillRegion( RegionID::Boundary_top );
    fillRegion( RegionID::Boundary_left );
    fillRegion( RegionID::Boundary_right );
}

// MESH::region::iterator
void structured2dRegularRectangle::fillRegion(RegionID id) {
    int nbX = static_cast<int>(nbCellsX_);
    int nbY = static_cast<int>(nbCellsY_);
    int nbXY = static_cast<int>(nbCells_);

    if (id == RegionID::Entire ) {
        regions[id] = Region{
            0,       // start
            1,       // step
            nbXY     // count
        };
    }else if ( id == RegionID::Boundary_left )
        regions[id] = Region{0,nbX,nbY  };
    else if ( id == RegionID::Boundary_bottom )
        regions[id] = Region{nbX*(nbY-1),1,nbX  };
    else if ( id == RegionID::Boundary_right )
        regions[id] = Region{nbX-1 +(nbY-1)*nbX,-nbX,nbY  };
    else if ( id == RegionID::Boundary_top )
        regions[id] = Region{nbX-1,-1,nbX  };
}

bool structured2dRegularRectangle::isBoundaryCell( unsigned int cell ) const {
        unsigned int i = cell % nbCellsX();
        unsigned int j = cell / nbCellsX();
        return i == 0 || i == nbCellsX() - 1 || j == 0 || j == nbCellsY() - 1;
}
const GLOBAL::scalar structured2dRegularRectangle::getCellSpacing_X( ) const
{
    return lenX_/static_cast<GLOBAL::scalar>(nbCellsX_);
}
const GLOBAL::scalar structured2dRegularRectangle::getCellSpacing_Y( ) const
{
    return lenY_/static_cast<GLOBAL::scalar>(nbCellsY_);
}

const GLOBAL::scalar structured2dRegularRectangle::getCellReciprocalSpacing_X( ) const
{
    return static_cast<GLOBAL::scalar>(nbCellsX_)/lenX_;
}

const GLOBAL::scalar structured2dRegularRectangle::getCellReciprocalSpacing_Y( ) const
{
    return static_cast<GLOBAL::scalar>(nbCellsY_)/lenY_;
}

const GLOBAL::scalar structured2dRegularRectangle::getCellFaceAreal_Y( ) const
{
    return getCellSpacing_Y()*cellThickness;
}

const GLOBAL::scalar structured2dRegularRectangle::getCellFaceAreal_X( ) const
{
    return getCellSpacing_X()*cellThickness;
}

const GLOBAL::scalar structured2dRegularRectangle::getCellCenterCoordinate_X(int cellId) const{
    unsigned int i = cellId % nbCellsX_;
    // Cell-centered coordinates: offset by +0.5 cell widths
    return (GLOBAL::scalar( i) + GLOBAL::scalar(0.5)) * getCellSpacing_X();;
}

const GLOBAL::scalar structured2dRegularRectangle::getCellCenterCoordinate_Y(int cellId) const{
    int i = cellId / nbCellsX_;
    return (nbCellsY()-i -0.5) * getCellSpacing_Y() ;
}

const sCoordinates structured2dRegularRectangle::getCellFacePos(RegionID id, int cellID) const
{
    switch (id)
    {
        case RegionID::Boundary_bottom:
            return sCoordinates(getCellCenterCoordinate_X(cellID),getCellCenterCoordinate_Y(cellID)-getCellSpacing_Y()*0.5);
        case RegionID::Boundary_top:
            return sCoordinates(getCellCenterCoordinate_X(cellID),getCellCenterCoordinate_Y(cellID)+getCellSpacing_Y()*0.5);
        case RegionID::Boundary_left:
            return sCoordinates(getCellCenterCoordinate_X(cellID)-getCellSpacing_X()*0.5,getCellCenterCoordinate_Y(cellID));
        case RegionID::Boundary_right:
            return sCoordinates(getCellCenterCoordinate_X(cellID)+getCellSpacing_X()*0.5,getCellCenterCoordinate_Y(cellID));
        default :
            return sCoordinates(0.0,0.0);
    }
}

