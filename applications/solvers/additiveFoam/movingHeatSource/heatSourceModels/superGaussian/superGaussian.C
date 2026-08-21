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

#include "superGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(superGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, superGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::superGaussian::superGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    s_(vector::zero)
{
    k_ = heatSourceModelCoeffs_.lookup<scalar>("k");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModels::superGaussian::update()
{
    updateDimensions();

    s_ = dimensions_/Foam::pow(2.0, 1.0/k_);

    const scalar fMax =
        Foam::pow
        (
            invIncGammaRatio_P(3.0/k_, 1.0 - profileTol_),
            1.0/k_
        );

    sourceMin_ = vector(-fMax*s_.x(), -fMax*s_.y(), -fMax*s_.z());

    sourceMax_ = vector(fMax*s_.x(), fMax*s_.y(), 0.0);

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            (2.0/3.0)*s_.x()*s_.y()*s_.z()*pi
           *Foam::tgamma(1.0 + 3.0/k_)
        );
}


inline Foam::scalar
Foam::heatSourceModels::superGaussian::weight(const vector& r) const
{
    if (r.z() > 0)
    {
        return 0;
    }

    return Foam::exp(-Foam::pow(magSqr(cmptDivide(r, s_)), k_/2.0));
}

bool Foam::heatSourceModels::superGaussian::read()
{
    if (heatSourceModel::read())
    {
        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("k") >> k_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
