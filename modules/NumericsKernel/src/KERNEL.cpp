#include "../include/KERNEL.h"
#include "LinEqsSolvers.h"
#include "blaze/Blaze.h"
#include "algorithm"
#include <ranges>

// make object registry to have own files.
KERNEL::ObjectRegistry::ObjectRegistry() = default;
KERNEL::ObjectRegistry::~ObjectRegistry() = default;
KERNEL::ObjectRegistry::ObjectRegistry(ObjectRegistry&&) noexcept = default;
KERNEL::ObjectRegistry& KERNEL::ObjectRegistry::operator=(ObjectRegistry&&) noexcept = default;

template<typename MT>
bool checkMatrixTypeIsSparse(const MT& A)
{
    if constexpr (blaze::IsSparseMatrix_v<MT>)
    {
        return true;
    }
    return false;
}

std::vector<int> KERNEL::getBandIDs(const KERNEL::smatrix &A){
    std::vector<int> bandIDs;

    // collecting indices of bands in first row:
    for (auto it = A.cbegin(0); it != A.cend(0); ++it) {
        bandIDs.push_back(static_cast<int>( it->index() ) );
    }

    // inserting (-)rowindex, if first nonzero entry is in first column.
    for( size_t r=1; r<A.rows(); ++r ) {
        auto it=A.cbegin(r);
        if ( it->index()==0 ) {
            bandIDs.insert(bandIDs.begin(), -static_cast<int>( r ) );
        }
    }

    return bandIDs;
}

std::vector<size_t> calcNonzeros(std::size_t N, const std::vector<int>& bandIDs) {

    // nonzeros are filled with number of bands in every row.
    std::vector<size_t> nonzeros(N, bandIDs.size());

    // For a negative band id, M, all rows [ 0, |M| [ are decremented.
    // For a positive band id, M, all rows [N-ri, M [ are decremented.
    for (auto bid : bandIDs) {
        if (bid<0) {
            for (auto& n : nonzeros | std::views::take(abs(bid))) --n;
        }else {
            for (auto& n : nonzeros | std::views::drop(N-bid)) --n;
        }
    }

    return nonzeros;
}


KERNEL::smatrix KERNEL::newTempBandedSMatrix(std::size_t N, const std::vector<int>& bandIDs, GLOBAL::scalar init) {

    auto entriesUnique = [](std::vector<int> v) {
        std::ranges::sort(v);
        return std::ranges::adjacent_find(v) == v.end();
    };

    if (!entriesUnique(bandIDs))
        throw std::runtime_error("bandIDs contains duplicates");

    // a vector of number of nonzero elements per matrix row
    std::vector<size_t> nonzeros = calcNonzeros(N, bandIDs);

    smatrix A(N,N,nonzeros);

    for (auto bID: bandIDs) {
        fillBand(blaze::band(A,bID), init);
    }

    return A;
}

// Vector creation
KERNEL::VectorHandle KERNEL::ObjectRegistry::newVector(size_t size, GLOBAL::scalar initialValue) {
    if (registryClosed_)
        throw std::runtime_error("Registry closed. New objects must be defined before closing registry.");

    auto id = nextID_++;

    registry_[id] = std::make_unique<KERNEL::vector>(size, initialValue);

    return VectorHandle{id};
}

// FJA combine this one with newBandedTemp function
KERNEL::MatrixHandle KERNEL::ObjectRegistry::newMatrix(size_t rows, size_t cols, bool sparse) {
    if (registryClosed_)
        throw std::runtime_error("Registry closed. New objects must be defined before closing registry.");

    auto id = nextID_++;
    // nextID_ = nextID_ + 1;
    // auto id = nextID_;

    if (sparse)
        registry_[id] = std::make_unique<KERNEL::smatrix>(rows, cols);
    else
        registry_[id] = std::make_unique<KERNEL::dmatrix>(rows, cols);

    return MatrixHandle{id};
}

KERNEL::vector& KERNEL::ObjectRegistry::getVectorRef(VectorHandle handle) {
    if (!registryClosed_) {
        throw std::runtime_error("Close registry before accessing objects.");
    }
    auto it = registry_.find(handle.id);
    if (it == registry_.end()) {
        throw std::runtime_error("Invalid object ID");
    }

    // Try to get vector from variant
    std::unique_ptr<KERNEL::vector>* vecPtr = std::get_if< std::unique_ptr<KERNEL::vector>>(&it->second);
    if (!vecPtr) {
        throw std::runtime_error("Object is not a vector");
    }

    return **vecPtr;
}

KERNEL::dmatrix& KERNEL::ObjectRegistry::getDenseMatrixRef(MatrixHandle handle) {

    auto it = registry_.find(handle.id);

    if (it == registry_.end())
        throw std::runtime_error("Invalid object ID");

    auto* matPtr = std::get_if<std::unique_ptr<KERNEL::dmatrix>>(&it->second);
    if (!matPtr)
        throw std::runtime_error("Object is not a dense matrix");

    return **matPtr;
}

KERNEL::smatrix& KERNEL::ObjectRegistry::getSparseMatrixRef(MatrixHandle handle) {

    auto it = registry_.find(handle.id);

    if (it == registry_.end())
        throw std::runtime_error("Invalid object ID");

    auto* matPtr = std::get_if<std::unique_ptr<KERNEL::smatrix>>(&it->second);
    if (!matPtr)
        throw std::runtime_error("Object is not a sparse matrix");

    return **matPtr;
}

// FJA assert that input matrix is dmatrix type, if not abort
void KERNEL::solve(const dmatrix& A, vector& x, const vector& b, KERNEL::SolverMethod method, const GLOBAL::scalar tolerance, const unsigned int maxIter) {

    if (method == BiCGSTAB) {
        LINEQSOLVERS::solve_BiCGSTAB(A, x, b, tolerance, maxIter);
    }else if (method == GaussSeidel) {
        LINEQSOLVERS::solve_GaussSeidel(A, x, b, tolerance, maxIter);
    }else if (method == Jacobi) {
        LINEQSOLVERS::solve_Jacobi(A, x, b, tolerance, maxIter);
    }else if (method == Blaze_automatic) {
        blaze::solve(A, x, b);
    }
}

void KERNEL::solve(const KERNEL::smatrix& A, KERNEL::vector& x, const KERNEL::vector& b,  KERNEL::SolverMethod method, const GLOBAL::scalar tolerance, const unsigned int maxIter) {

    // not here!
    // checkLinEqSystemConsistency(A,b);

    if (method == BiCGSTAB) {
        LINEQSOLVERS::solve_BiCGSTAB(A, x, b, tolerance, maxIter);
    }else {
        throw std::runtime_error("WARNING: the sparse matrix cannot be solved with the chosen linear solver.");
    }

}
