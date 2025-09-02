
#include "denseMatrixSolvers.h"

#include <iostream>
#include "MathUtils.h"

std::vector<double> Dense::GaussSeidel::solve()
{
    return std::vector<double>(b_.size());
}

std::vector<double> Dense::JacobiIter::solve()
{
    std::vector<double> x, x_old;
    const auto n = b_.size();
    x_old.resize(n);
    x.resize(n);
    x_old.assign(n, 0);
    double sum_terms = 0.0;
    for ( int k = 0; k < maxIter_; k++)
    {
        for (int i = 0; i < n; i++)
        {
            sum_terms = 0.0;

            for (int j = 0; j < i; j++)
            {
                sum_terms += A_[i][j] * x_old[j];
            }

            for (int j = i + 1; j < n; j++)
            {
                sum_terms += A_[i][j] * x_old[j];
            }

            x[i] = (-sum_terms + b_[i] ) / A_[i][i];
        }

        if ( MathUtils::norm2(x,x_old) < tolerance_)
        {
            std::cout<<"Jacobi solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
            return x;
       }
        x_old = x;
    }
    std::cerr<<"Jacobi solver did not converge within"<<  std::to_string(maxIter_)<<" iterations."<<std::endl;
    return x;
}

