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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace heatSourceProjections
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable
    (
        heatSourceProjection,
        exponential,
        dictionary
    );
}
}

Foam::heatSourceProjections::exponential::exponential
(
    const dictionary& dict
)
:
    nSlope_(0),
    nIntercept_(0),
    d_(0),
    k_(1),
    integral_(0),
    zMin_(0)
{
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");

    nSlope_ = coeffs.lookup<scalar>("nSlope");
    nIntercept_ = coeffs.lookup<scalar>("nIntercept");
}


void Foam::heatSourceProjections::exponential::update
(
    const scalar depth,
    const scalar aspectRatio,
    const scalar tolerance
)
{
    d_ = depth;

    const scalar n =
        min
        (
            max
            (
                nSlope_*std::log2(aspectRatio) + nIntercept_,
                0.0
            ),
            9.0
        );

    k_ = Foam::pow(2.0, n);

    integral_ =
        d_*Foam::tgamma(1.0/k_)
       /(k_*Foam::pow(3.0, 1.0/k_));

    zMin_ =
       -d_
       *Foam::pow
        (
            invIncGammaRatio_P(1.0/k_, 1.0 - tolerance)/3.0,
            1.0/k_
        );
}


inline Foam::scalar Foam::heatSourceProjections::exponential::weight
(
    const scalar z
) const
{
    if (z > 0)
    {
        return 0;
    }

    return Foam::exp(-3.0*Foam::pow(-z/d_, k_));
}

// ************************************************************************* //
