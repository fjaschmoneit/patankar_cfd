//
// Created by Peter Berg Ammundsen on 02/09/2025.
//

#pragma once
#include <vector>

class MathUtils
{
public:
    // ||x - y||_2  (Euclidean norm of the difference)
    static double diffNorm2(const std::vector<double>& x,
                        const std::vector<double>& y);
private:
    static double norm2(const std::vector<double>& x);
};


