//
// Created by Peter Berg Ammundsen on 30/09/2025.
//

#ifndef FVMCalculatedCoefficient_H
#define FVMCalculatedCoefficient_H
#include <globalTypeDefs.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <ranges>


struct FvmCoeffs
{
    // Geometri
    unsigned int Nx_, Ny_;
    unsigned int numberOfCells_;
    GLOBAL::scalar dx_, dy_;
    GLOBAL::scalar L_;

    // Coefficients (per cell)
    GLOBAL::vector aP_, aW_, aE_, aS_, aN_, Su_, Sp_;

    // Diffusion conductances (D) and distances (d) per direction
    //GLOBAL::scalar Dw_, De_, Ds_, Dn_;
    GLOBAL::scalar dW_, dE_, dS_, dN_;

    GLOBAL::scalar k_;
    GLOBAL::scalar A_;
   
};

class FVMCalculateCoefficient
{
public:
    FvmCoeffs fvmCoeffs{};
    FVMCalculateCoefficient(unsigned int Nx,unsigned int Ny,unsigned int numberOfCells,GLOBAL::scalar Length)
    {
        fvmCoeffs.Nx_ = Nx;
        fvmCoeffs.Ny_ = Ny;
        fvmCoeffs.numberOfCells_ = numberOfCells;
        fvmCoeffs.L_ = Length;
        fvmCoeffs.dx_ = Length/Nx;
        fvmCoeffs.dy_ = Length/Ny;
        fvmCoeffs.k_ = 200;
        fvmCoeffs.A_ = 10e-3;
        fvmCoeffs.aE_.resize(numberOfCells);
        fvmCoeffs.aS_.resize(numberOfCells);
        fvmCoeffs.aN_.resize(numberOfCells);
        fvmCoeffs.aW_.resize(numberOfCells);
        fvmCoeffs.aP_.resize(numberOfCells);
        fvmCoeffs.Su_.resize(numberOfCells);
        fvmCoeffs.Sp_.resize(numberOfCells);

        fvmCoeffs.dW_ = 0;
        fvmCoeffs.dE_ = 0;

        fvmCoeffs.dS_ = 0;
        fvmCoeffs.dN_ = 0;
    };
    void saveAsCSV(const std::vector<double>& vec, const std::string& filename);
    void calculateDiffusionAidx(const FVM::CardinalDirection& dir)
    {
        auto value = fvmCoeffs.k_*fvmCoeffs.A_/fvmCoeffs.dx_;
        if (dir != FVM::CardinalDirection::centre)
        {
            addCoefficientsToWestEastNortOrSouth(dir, value);
        }
    }
    void calculateDiffusionAidy(const FVM::CardinalDirection& dir)
    {
        auto value = fvmCoeffs.k_*fvmCoeffs.A_/fvmCoeffs.dy_;
        if (dir != FVM::CardinalDirection::centre)
        {
            addCoefficientsToWestEastNortOrSouth(dir, value);
        }
    }

    [[nodiscard]] GLOBAL::vector getAi(FVM::CardinalDirection dir) const;
    [[nodiscard]] GLOBAL::vector getB() const;
    
    [[nodiscard]] std::vector<unsigned int> getFaceNumberSouth() const
    {
        std::vector<unsigned int> values(fvmCoeffs.Ny_);
        std::iota(values.begin(), values.end(), static_cast<int>(2 * fvmCoeffs.Nx_));
        return values;
    }
    [[nodiscard]] std::vector<unsigned int> getFaceNumberNorth() const
    {
        std::vector<unsigned int> values(fvmCoeffs.Ny_);
        std::iota(values.begin(), values.end(), static_cast<int>(0));
        return values;
    }
    [[nodiscard]] std::vector<unsigned int> getFaceNumberWest() const
    {
        std::vector<unsigned int> values(fvmCoeffs.Ny_);
        std::iota(values.begin(), values.end(), static_cast<int>(0));
        for (auto &x : values | std::views::all)
        {
            x *= fvmCoeffs.Nx_;
        }
        return values;
    }

    [[nodiscard]] std::vector<unsigned int> getFaceNumberEast() const
    {
        std::vector<unsigned int> values(fvmCoeffs.Ny_);
        std::iota(values.begin(), values.end(), static_cast<int>(0));
        for (auto &x : values | std::views::all)
        {
            x = x * (fvmCoeffs.Nx_)+ fvmCoeffs.Ny_-1;
        }
        return values;
    }

    void apply_Boundary(const FVM::CardinalDirection& dir,const GLOBAL::vector& value)
    {
        std::vector<unsigned int> Faces;
        GLOBAL::vector* ai = &fvmCoeffs.aE_;
        switch (dir)
        {
            case FVM::CardinalDirection::east:  Faces = getFaceNumberEast(); ai = &fvmCoeffs.aE_; break;
            case FVM::CardinalDirection::south: Faces = getFaceNumberSouth();ai = &fvmCoeffs.aS_; break;
            case FVM::CardinalDirection::north: Faces = getFaceNumberNorth();ai = &fvmCoeffs.aN_; break;
            case FVM::CardinalDirection::west:  Faces = getFaceNumberWest(); ai = &fvmCoeffs.aW_; break;
            default: throw std::invalid_argument("Unknown direction");
        }

        unsigned int i = 0;
        //apply Direchtly bounadry condition.
        for (const unsigned int face : Faces)
        {
            fvmCoeffs.Sp_[face] += -2 * (*ai)[face];
            fvmCoeffs.Su_[face] +=  2 * (*ai)[face]*value[i];
            (*ai)[face] = 0;
            i++;
        }
    }

    void CollectAp()
    {
        const auto N = static_cast<std::size_t>(fvmCoeffs.numberOfCells_);

        for (auto p = 0; p < N; ++p)
        {
            fvmCoeffs.aP_[p] =
                  fvmCoeffs.aE_[p]
                + fvmCoeffs.aW_[p]
                + fvmCoeffs.aS_[p]
                + fvmCoeffs.aN_[p]
                - fvmCoeffs.Sp_[p];
        }
    }


private:

    void addCoefficientsToWestEastNortOrSouth(const FVM::CardinalDirection& dir, GLOBAL::scalar& ai);
};



#endif //FVMCalculatedCoefficient_H
