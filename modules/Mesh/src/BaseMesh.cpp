//
// Created by Peter Berg Ammundsen on 21/01/2026.
//

#include "BaseMesh.h"
using namespace MESH;
BaseMesh::BaseMesh() {
}

const Region& BaseMesh::region(RegionID id) const
{
    auto it = regions.find(id);
    if (it == regions.end()) throw std::out_of_range("Unknown region id");
    return it->second;
}