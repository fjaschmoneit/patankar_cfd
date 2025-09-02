//
// Created by Peter Berg Ammundsen on 02/09/2025.
//

#include "MathUtils.h"
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>


double MathUtils::norm2(const std::vector<double>& x,
                    const std::vector<double>& y)
{
    auto dq = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    return std::sqrt(dq);
}