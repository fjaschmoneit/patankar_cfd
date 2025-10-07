//
// Created by Peter Berg Ammundsen on 30/09/2025.
//

#include "FVMCalculateCoefficient.h"
#include <functional>
#include <algorithm>
#include <fstream>

GLOBAL::vector FVMCalculateCoefficient::getAi(FVM::CardinalDirection dir) const
{
    switch (dir) {
        case FVM::CardinalDirection::centre: return fvmCoeffs.aP_;
        case FVM::CardinalDirection::east:   return fvmCoeffs.aE_;
        case FVM::CardinalDirection::south:  return fvmCoeffs.aS_;
        case FVM::CardinalDirection::west:   return fvmCoeffs.aW_;
        case FVM::CardinalDirection::north:  return fvmCoeffs.aN_;
        default: throw std::invalid_argument("Unknown direction");
    }
}
GLOBAL::vector FVMCalculateCoefficient::getB() const
{
    return fvmCoeffs.Su_;
}

void FVMCalculateCoefficient::addCoefficientsToWestEastNortOrSouth(const FVM::CardinalDirection& dir, GLOBAL::scalar& value)
{
    switch (dir)
    {
        case FVM::CardinalDirection::east:   std::fill(fvmCoeffs.aE_.begin(), fvmCoeffs.aE_.end(), value); break;
        case FVM::CardinalDirection::south:  std::fill(fvmCoeffs.aS_.begin(), fvmCoeffs.aS_.end(), value); break;
        case FVM::CardinalDirection::west:   std::fill(fvmCoeffs.aW_.begin(), fvmCoeffs.aW_.end(), value); break;
        case FVM::CardinalDirection::north:  std::fill(fvmCoeffs.aN_.begin(), fvmCoeffs.aN_.end(), value); break;
        default: throw std::invalid_argument("Unknown direction");
    }
}

void FVMCalculateCoefficient::saveAsCSV(const std::vector<double>& vec, const std::string& filename) {
    std::ofstream file(filename);
    for (size_t i = 0; i < vec.size(); ++i)
    {
        file << vec[i] << "\n";
    }
    file.close();
    //saveAsCSV(field,"/Users/ri03jm/Library/Mobile Documents/com~apple~CloudDocs/DropboxFolder/Peter/skole/PhD/MatlabScript/output.csv");
}
