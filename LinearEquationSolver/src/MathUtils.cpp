//
// Created by Peter Berg Ammundsen on 02/09/2025.
//

#include "MathUtils.h"
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>


double MathUtils::norm2(const std::vector<double>& x)
{
    auto dq = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    return std::sqrt(dq);
}

double MathUtils::diffNorm2(const std::vector<double>& x,const std::vector<double>& y)
{
    std::vector<double> diff(x.size());
    std::transform(x.begin(), x.end(), y.begin(), diff.begin(),
                   [](double a, double b){ return a - b; });
    return norm2(diff);

}