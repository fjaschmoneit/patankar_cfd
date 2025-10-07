//
// Created by Peter Berg Ammundsen on 02/09/2025.
//

#pragma once
#include <vector>
#include <blaze/Math.h>

class MathUtils
{
public:
    // ||x - y||_2  (Euclidean norm of the difference)
    template <typename VectorType>
    static double diffNorm2(const VectorType& x, const VectorType& y)
    {
        VectorType diff(x.size());
        std::transform(x.begin(), x.end(), y.begin(), diff.begin(),
                       [](double a, double b){ return a - b; });
        return norm2(diff);
    }


private:
    template <typename VectorType>
    static double norm2(const VectorType& x)
    {
        auto dq = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
        return std::sqrt(dq);
    }
};


