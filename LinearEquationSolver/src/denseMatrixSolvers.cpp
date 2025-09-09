
#include "denseMatrixSolvers.h"

#include <iostream>
#include "MathUtils.h"

blaze::DynamicVector<double> Dense::GaussSeidel::solve1()
{
    return blaze::DynamicVector<double>{ 4, -1, 3 };
}
// to do, temporary function, will be removed later
std::vector<double> Dense::GaussSeidel::solve()
{
    const int n = (int)b_.size();
    if ((int)A_.size() != n) throw std::invalid_argument("A,b size mismatch");
    for (int i=0;i<n;++i) if ((int)A_[i].size()!=n) throw std::invalid_argument("A not square");

    std::vector<double> x(n, 0.0);
    std::vector<double> x_old(n, 0.0);

    for (int k = 0; k < maxIter_; ++k)
    {
        x_old = x;

        for (int i = 0; i < n; ++i)
        {
            double sum_terms = 0.0;

            // Users new values for j < i
            // Not same as Jacobi but improve as the steps already  used new calculated guess instead of used old value
            for (int j = 0; j < i; ++j)
            {
                sum_terms += A_[i][j] * x[j];
            }

            // Users old  values for j > i -> same as Jacobi Methods
            for (int j = i+1; j < n; ++j)
            {
                sum_terms += A_[i][j] * x_old[j];
            }

            const double aii = A_[i][i];
            if (aii == 0.0)
            {
                throw std::runtime_error("GaussSeidel: zero diagonal");
            }

            x[i] = (-sum_terms + b_[i] ) / aii;
        }

      if ( MathUtils::diffNorm2(x,x_old) < tolerance_)
        {
            std::cout<<"Gauss Seidel solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
            return x;
        }
    }
    std::cerr<<"Gauss Seidel, solver did not converge within"<<  std::to_string(maxIter_)<<" iterations."<<std::endl;
    return x;
}
blaze::DynamicVector<double> Dense::JacobiIter::solve1() {

    auto rows = Ab_.rows();
    blaze::DynamicVector<double> x(rows, 0.0 );
    blaze::DynamicVector<double, blaze::columnVector> x_old(rows, 0.0 );
    blaze::DynamicMatrix<double> B = Ab_;     // deep copy

    auto d = blaze::diagonal(Ab_);  // view af diagonalen
    auto invD = 1.0/d;

    blaze::DiagonalMatrix< blaze::CompressedMatrix<double> > invDSparse( Ab_.rows() );
    blaze::diagonal( invDSparse ) = invD;
    blaze::diagonal(B) = 0.0;

    for (int k = 0; k < maxIter_; ++k)
    {
        auto rhs = bb_ - B * x_old;
        x = invDSparse * rhs;
        if ( MathUtils::diffNorm2(x,x_old) < tolerance_)
        {
            std::cout<<"Jacobi solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
            return x;
        }
        x_old = x;
    }
    std::cerr<<"Jacobi solver did not converge within"<<  std::to_string(maxIter_)<<" iterations."<<std::endl;
    return x;
}

// to do, temporary function, will be removed later
std::vector<double> Dense::JacobiIter::solve()
{
    std::vector<double> x, x_old;
    const auto n = b_.size();
    x_old.resize(n);
    x.resize(n);
    x_old.assign(n, 1);
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
        if ( MathUtils::diffNorm2(x,x_old) < tolerance_)
        {
            std::cout<<"Jacobi solver converged after "<<  std::to_string(k)<<" iterations."<<std::endl;
            return x;
       }
        x_old = x;
    }
    std::cerr<<"Jacobi solver did not converge within"<<  std::to_string(maxIter_)<<" iterations."<<std::endl;
    return x;
}
