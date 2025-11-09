#include "../include/KERNEL.h"
// #include "LinEqsSolvers.h"
#include "blaze/Blaze.h"

ObjectRegistry::ObjectRegistry() = default;
ObjectRegistry::~ObjectRegistry() = default;
ObjectRegistry::ObjectRegistry(ObjectRegistry&&) noexcept = default;
ObjectRegistry& ObjectRegistry::operator=(ObjectRegistry&&) noexcept = default;


// Vector creation
VectorHandle ObjectRegistry::newVector(size_t size, KERNEL::scalar initialValue) {
    if (registryClosed_)
        throw std::runtime_error("Registry closed. New objects must be defined before closing registry.");

    auto id = nextID_++;

    registry_[id] = std::make_unique<KERNEL::vector>(size, initialValue);

    return VectorHandle{id};
}

MatrixHandle ObjectRegistry::newMatrix(size_t rows, size_t cols, bool sparse) {
    if (registryClosed_)
        throw std::runtime_error("Registry closed. New objects must be defined before closing registry.");

    auto id = nextID_++;

    if (sparse)
        registry_[id] = std::make_unique<KERNEL::smatrix>(rows, cols);
    else
        registry_[id] = std::make_unique<KERNEL::dmatrix>(rows, cols);

    return MatrixHandle{id};
}

KERNEL::vector& ObjectRegistry::getVectorRef(VectorHandle handle) {
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

KERNEL::dmatrix& ObjectRegistry::getDenseMatrixRef(MatrixHandle handle) {

    auto it = registry_.find(handle.id);

    if (it == registry_.end())
        throw std::runtime_error("Invalid object ID");

    auto* matPtr = std::get_if<std::unique_ptr<KERNEL::dmatrix>>(&it->second);
    if (!matPtr)
        throw std::runtime_error("Object is not a dense matrix");

    return **matPtr;
}

KERNEL::smatrix& ObjectRegistry::getSparseMatrixRef(MatrixHandle handle) {

    auto it = registry_.find(handle.id);

    if (it == registry_.end())
        throw std::runtime_error("Invalid object ID");

    auto* matPtr = std::get_if<std::unique_ptr<KERNEL::smatrix>>(&it->second);
    if (!matPtr)
        throw std::runtime_error("Object is not a sparse matrix");

    return **matPtr;
}