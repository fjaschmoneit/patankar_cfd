#include <vector>
#include "modules/NumericsKernel/include/KERNEL.h"

std::vector<int> bandIDs{-3,0,3,5};

auto A = KERNEL::newTempBandedSMatrix(7, bandIDs);

