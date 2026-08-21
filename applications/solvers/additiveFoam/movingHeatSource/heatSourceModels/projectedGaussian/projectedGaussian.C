/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023-2026 Oak Ridge National Laboratory
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "projectedGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(projectedGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, projectedGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::projectedGaussian::projectedGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    k_(scalar(1))
{
    nSlope_ = heatSourceModelCoeffs_.lookup<scalar>("nSlope");

    nIntercept_ =
        heatSourceModelCoeffs_.lookup<scalar>("nIntercept");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModels::projectedGaussian::update()
{
    updateDimensions();

    const scalar a =
        dimensions_.z()
       /min(staticDimensions_.x(), staticDimensions_.y());

    const scalar n =
        min(max(nSlope_*std::log2(a) + nIntercept_, 0.0), 9.0);

    k_ = std::pow(2.0, n);

    // (1 - epsilon)^2 = (1 - profileTol_)
    const scalar epsilon =
        profileTol_/(1.0 + Foam::sqrt(1.0 - profileTol_));

    const scalar fMax = Foam::sqrt(-0.5*std::log(epsilon));

    const scalar zMax =
        dimensions_.z()
       *Foam::pow
        (
            invIncGammaRatio_P(1.0/k_, 1.0 - epsilon)/3.0,
            1.0/k_
        );

    sourceMin_ =
        vector(-fMax*dimensions_.x(), -fMax*dimensions_.y(), -zMax);

    sourceMax_ =
        vector(fMax*dimensions_.x(), fMax*dimensions_.y(), 0);

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            0.5*pi*dimensions_.x()*dimensions_.y()*dimensions_.z()
           *Foam::tgamma(1.0/k_)
           /(k_*std::pow(3.0, 1.0/k_))
        );
}

inline Foam::scalar
Foam::heatSourceModels::projectedGaussian::weight(const vector& r) const
{
    if (r.z() > 0)
    {
        return 0;
    }

    const scalar I =
        std::exp
        (
            -2.0
          * (
                Foam::sqr(r.x()/dimensions_.x())
              + Foam::sqr(r.y()/dimensions_.y())
            )
        );

    return I*std::exp(-3.0*std::pow(-r.z()/dimensions_.z(), k_));
}

bool Foam::heatSourceModels::projectedGaussian::read()
{
    if (heatSourceModel::read())
    {
        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("nSlope") >> nSlope_;
        heatSourceModelCoeffs_.lookup("nIntercept")
            >> nIntercept_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
