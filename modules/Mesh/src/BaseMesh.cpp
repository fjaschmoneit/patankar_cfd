//
// Created by Peter Berg Ammundsen on 21/01/2026.
//

#include "BaseMesh.h"
using namespace MESH;
BaseMesh::BaseMesh(GLOBAL::scalar lengthX, unsigned int nbCellsX, GLOBAL::scalar lengthY, unsigned int nbCellsY)
    : lenX_(lengthX)
    , lenY_(lengthY)
    , nbCellsX_(nbCellsX)
    , nbCellsY_(nbCellsY)
    , nbCells_(nbCellsX * nbCellsY)
{}
unsigned int BaseMesh::nbCellsX() const
{
    return nbCellsX_;
}
unsigned int BaseMesh::nbCellsY() const
{
    return nbCellsY_;
}
unsigned int BaseMesh::nbCells()  const
{
    return nbCells_;
}
GLOBAL::scalar BaseMesh::lenX() const
{
    return lenX_;
}
GLOBAL::scalar BaseMesh::lenY() const
{
    return lenY_;
}

const Region& BaseMesh::region(RegionID id) const
{
    auto it = regions.find(id);
    if (it == regions.end()) throw std::out_of_range("Unknown region id");
    return it->second;
}