//
// Created by Peter Berg Ammundsen on 08/12/2025.
//

#ifndef STRUCTURED2D_H
#define STRUCTURED2D_H
#pragma once
#include <map>
#include "GlobalTypeDefs.h"
#include "BaseMesh.h"


namespace MESH {

    class structured2dRegularRectangle: public BaseMesh
    {
        // Registry of all regions; enum class ensures type safety
        //std::unordered_map<RegionID, Region> regions;

        // how can I do this at compile time?
        // Each boundary side is represented by a vector of cell indices
        // example: rectangular mesh with (nbX,nbY)=(5,4) cells
        // its cell and boundary indices are illustrated here:
        //   y ^
        //     |          3
        //     +----+----+----+----+----+
        //     |  0 |  1 |  2 |  3 |  4 |
        //     +----+----+----+----+----+
        //     |  5 |  6 |  7 |  8 |  9 |
        //  0  +----+----+----+----+----+    2
        //     | 10 | 11 | 12 | 13 | 14 |
        //     +----+----+----+----+----+
        //     | 15 | 16 | 17 | 18 | 19 |
        //     +----+----+----+----+----+ --> x
        //               1
        //
        // the boundaries are names counter-clockwise, starting with the left one.
        // the cell ids per boundary also follow a counter-clockwise direction.
        // boundaries[0] = [0,5,10,15]
        // boundaries[1] = [15,16,17,18,19]
        // boundaries[2] = [19,14,9,4]
        // boundaries[3] = [4,3,2,1,0]
    private:
        void fillRegion(RegionID region) override;
        // order must follow initialization list
        const GLOBAL::scalar lenX_;
        const GLOBAL::scalar lenY_;
        const unsigned int nbCellsX_;
        const unsigned int nbCellsY_;
        const unsigned int nbCells_;

    public:
        structured2dRegularRectangle(GLOBAL::scalar lengthX, unsigned int nbCellsX, GLOBAL::scalar lengthY, unsigned int nbCellsY);
        void init();
        unsigned int nbCellsX() const;
        unsigned int nbCellsY() const;
        unsigned int nbCells()  const;
        GLOBAL::scalar lenX() const;
        GLOBAL::scalar lenY() const;

        const GLOBAL::scalar getCellSpacing_X( ) const;
        const GLOBAL::scalar getCellSpacing_Y( ) const;

        // returns a cell 1/dx and 1/dy.
        // in regular meshes, these are same everywhere
        // should accept iterator to mesh cell
        const GLOBAL::scalar getCellReciprocalSpacing_X( ) const;
        const GLOBAL::scalar getCellReciprocalSpacing_Y( ) const;
        const GLOBAL::scalar getCellFaceArea_Y( ) const;
        const GLOBAL::scalar getCellFaceArea_X( ) const;

        //virtual functions
        const GLOBAL::scalar getCellCenterCoordinate_X(int cellId) const override;
        const GLOBAL::scalar getCellCenterCoordinate_Y(int cellId) const override;

        bool isBoundaryCell( unsigned int i ) const override;
        const std::map<int,Coordinate> getCellFacesPos(RegionID id) const override;
    };

}
#endif //STRUCTURED2D_H
