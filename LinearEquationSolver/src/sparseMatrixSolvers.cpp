#include "sparseMatrixSolvers.h"

std::vector<double> Sparse::GaussSeidel::solve(){
    return std::vector<double>(b_.size(), 0.0);
}
