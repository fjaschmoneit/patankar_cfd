//
// Created by Peter Berg Ammundsen on 02/09/2025.
//

#include "MathUtils.h"
#include <cmath>
#include <stdexcept>
#include <vector>


double MathUtils::norm2(const std::vector<double>& x,
                    const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Math::Utils::norm2: vectors must have same size");
    }
    double sum_sq = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double d = x[i] - y[i];
        sum_sq += d * d;
    }
    return std::sqrt(sum_sq);
}