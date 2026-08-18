/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tabulatedProfile.H"
#include "IOstreams.H"
#include "error.H"

#include <cmath>
#include <fstream>
#include <limits>
#include <string>

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulatedProfile::tabulatedProfile()
:
    nx_(0),
    ny_(0),
    x0_(0),
    y0_(0),
    dx_(0),
    dy_(0),
    x1_(0),
    y1_(0),
    invDx_(0),
    invDy_(0),
    values_(0),
    integral_(0),
    centroidX_(0),
    centroidY_(0),
    d4SigmaMajor_(0),
    d4SigmaMinor_(0),
    d4SigmaEquivalent_(0),
    azimuth_(0)
{}

Foam::tabulatedProfile::tabulatedProfile(const fileName& profileFile)
:
    tabulatedProfile()
{
    read(profileFile);
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::tabulatedProfile::integrate()
{
    integral_ = 0;

    scalar firstMomentX = 0;
    scalar firstMomentY = 0;
    scalar secondMomentX = 0;
    scalar secondMomentY = 0;
    scalar crossMoment = 0;

    for (label j=0; j<ny_ - 1; ++j)
    {
        for (label i=0; i<nx_ - 1; ++i)
        {
            const label index = i + nx_*j;
            const scalar cellX = x0_ + i*dx_;
            const scalar cellY = y0_ + j*dy_;

            const scalar xWeights[3][2] =
            {
                {0.5*dx_, 0.5*dx_},
                {
                    dx_*(0.5*cellX + dx_/6.0),
                    dx_*(0.5*cellX + dx_/3.0)
                },
                {
                    dx_*
                    (
                        0.5*sqr(cellX) + cellX*dx_/3.0 + sqr(dx_)/12.0
                    ),
                    dx_*
                    (
                        0.5*sqr(cellX)
                      + 2.0*cellX*dx_/3.0
                      + sqr(dx_)/4.0
                    )
                }
            };

            const scalar yWeights[3][2] =
            {
                {0.5*dy_, 0.5*dy_},
                {
                    dy_*(0.5*cellY + dy_/6.0),
                    dy_*(0.5*cellY + dy_/3.0)
                },
                {
                    dy_*
                    (
                        0.5*sqr(cellY) + cellY*dy_/3.0 + sqr(dy_)/12.0
                    ),
                    dy_*
                    (
                        0.5*sqr(cellY)
                      + 2.0*cellY*dy_/3.0
                      + sqr(dy_)/4.0
                    )
                }
            };

            const scalar cellValues[2][2] =
            {
                {values_[index], values_[index + nx_]},
                {values_[index + 1], values_[index + nx_ + 1]}
            };

            for (label xNode=0; xNode<2; ++xNode)
            {
                for (label yNode=0; yNode<2; ++yNode)
                {
                    const scalar value = cellValues[xNode][yNode];

                    integral_ +=
                        value*xWeights[0][xNode]*yWeights[0][yNode];

                    firstMomentX +=
                        value*xWeights[1][xNode]*yWeights[0][yNode];

                    firstMomentY +=
                        value*xWeights[0][xNode]*yWeights[1][yNode];

                    secondMomentX +=
                        value*xWeights[2][xNode]*yWeights[0][yNode];

                    secondMomentY +=
                        value*xWeights[0][xNode]*yWeights[2][yNode];

                    crossMoment +=
                        value*xWeights[1][xNode]*yWeights[1][yNode];
                }
            }
        }
    }

    if (!std::isfinite(integral_) || integral_ <= 0)
    {
        FatalErrorInFunction
            << "Tabulated profile integral must be finite and positive"
            << exit(FatalError);
    }

    centroidX_ = firstMomentX/integral_;
    centroidY_ = firstMomentY/integral_;

    const scalar varianceX =
        secondMomentX/integral_ - sqr(centroidX_);

    const scalar varianceY =
        secondMomentY/integral_ - sqr(centroidY_);

    const scalar covariance =
        crossMoment/integral_ - centroidX_*centroidY_;

    if
    (
        !std::isfinite(varianceX)
     || !std::isfinite(varianceY)
     || !std::isfinite(covariance)
     || varianceX <= 0
     || varianceY <= 0
    )
    {
        FatalErrorInFunction
            << "Tabulated profile must have finite positive lateral width"
            << exit(FatalError);
    }

    const scalar meanVariance = 0.5*(varianceX + varianceY);
    const scalar varianceRadius =
        Foam::sqrt(0.25*sqr(varianceX - varianceY) + sqr(covariance));

    const scalar majorVariance = meanVariance + varianceRadius;
    const scalar minorVariance = meanVariance - varianceRadius;

    if (minorVariance <= 0)
    {
        FatalErrorInFunction
            << "Tabulated profile must have positive principal widths"
            << exit(FatalError);
    }

    d4SigmaMajor_ = 4.0*Foam::sqrt(majorVariance);
    d4SigmaMinor_ = 4.0*Foam::sqrt(minorVariance);
    d4SigmaEquivalent_ = Foam::sqrt(d4SigmaMajor_*d4SigmaMinor_);

    azimuth_ =
        varianceRadius
      <= 64.0*std::numeric_limits<scalar>::epsilon()*meanVariance
      ? 0.0
      : 0.5*Foam::atan2(2.0*covariance, varianceX - varianceY);
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tabulatedProfile::read(const fileName& profileFile)
{
    std::ifstream input(profileFile.c_str());

    if (!input.good())
    {
        FatalErrorInFunction
            << "Cannot open tabulated profile " << profileFile
            << exit(FatalError);
    }

    long long nx = 0;
    long long ny = 0;

    if (!(input >> nx >> ny) || nx < 2 || ny < 2)
    {
        FatalErrorInFunction
            << "Tabulated profile " << profileFile
            << " must begin with dimensions of at least 2 by 2"
            << exit(FatalError);
    }

    if
    (
        nx > std::numeric_limits<label>::max()
     || ny > std::numeric_limits<label>::max()
     || nx*ny > std::numeric_limits<label>::max()
    )
    {
        FatalErrorInFunction
            << "Tabulated profile " << profileFile << " is too large"
            << exit(FatalError);
    }

    nx_ = nx;
    ny_ = ny;

    if
    (
        !(input >> x0_ >> y0_ >> dx_ >> dy_)
     || !std::isfinite(x0_)
     || !std::isfinite(y0_)
     || !std::isfinite(dx_)
     || !std::isfinite(dy_)
     || dx_ <= 0
     || dy_ <= 0
    )
    {
        FatalErrorInFunction
            << "Tabulated profile " << profileFile
            << " requires finite coordinates and positive spacing"
            << exit(FatalError);
    }

    x1_ = x0_ + (nx_ - 1)*dx_;
    y1_ = y0_ + (ny_ - 1)*dy_;
    invDx_ = 1.0/dx_;
    invDy_ = 1.0/dy_;

    values_.setSize(nx_*ny_);

    forAll(values_, i)
    {
        if
        (
            !(input >> values_[i])
         || !std::isfinite(values_[i])
         || values_[i] < 0
        )
        {
            FatalErrorInFunction
                << "Tabulated profile " << profileFile
                << " requires exactly " << values_.size()
                << " finite nonnegative values"
                << exit(FatalError);
        }
    }

    std::string trailing;
    if (input >> trailing)
    {
        FatalErrorInFunction
            << "Unexpected trailing data in tabulated profile " << profileFile
            << exit(FatalError);
    }

    integrate();
}

// ************************************************************************* //
