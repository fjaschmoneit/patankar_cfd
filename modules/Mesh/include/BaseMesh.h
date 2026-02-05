//
// Created by Peter Berg Ammundsen on 21/01/2026.
//

#ifndef CBASEMESH_H
#define CBASEMESH_H
#include <cstddef>
#include <iterator>
#include <unordered_map>
#include "GlobalTypeDefs.h"
#include <map>

namespace MESH {
    // Enum for region identifiers
    enum class RegionID {
        Entire,
        Boundary_left,
        Boundary_bottom,
        Boundary_right,
        Boundary_top
    };

    struct sCoordinates
    {
        GLOBAL::scalar x = 0.0;
        GLOBAL::scalar y = 0.0;

        sCoordinates(GLOBAL::scalar _x, GLOBAL::scalar _y): x(_x), y(_y){}
    };

    struct Region {
        struct iterator {
            using value_type        = int;
            using difference_type   = std::ptrdiff_t;
            using reference         = int;            // yields by value
            using pointer           = const int*;     // not dereferenced normally
            using iterator_category = std::forward_iterator_tag;

            int current;
            int step;
            int remaining;

            reference operator*() const { return current; }
            iterator& operator++() { current += step; --remaining; return *this; }
            iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
            bool operator==(const iterator& other) const { return remaining == other.remaining; }
            bool operator!=(const iterator& other) const { return remaining != other.remaining; }
        };

        int start;
        int step;
        int count;

        iterator begin() const { return { start, step, count }; }
        iterator end()   const { return { 0, step, 0 }; }
        iterator last()  const { return { start + step*(count - 1), step, 1 }; }
    };

    class BaseMesh
    {
    protected:
        // Registry of all regions; enum class ensures type safety
        std::unordered_map<RegionID, Region> regions;

        // Derived meshes define how regions are populated
        virtual void fillRegion(RegionID region) = 0;

    public:
        BaseMesh();

        virtual ~BaseMesh() = default;

        const Region& region(RegionID id) const;
        GLOBAL::scalar cellThickness = 1;

        virtual const GLOBAL::scalar getCellCenterCoordinate_X(int cellId) const = 0;
        virtual const GLOBAL::scalar getCellCenterCoordinate_Y(int cellId) const = 0;

        virtual bool isBoundaryCell( unsigned int i ) const = 0;
        virtual const std::map<int,sCoordinates>  getCellFacesPos(RegionID id) const = 0;
    };
}

#endif //CBASEMESH_H
