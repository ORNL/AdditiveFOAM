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
#include "superGaussianProfile.H"
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
    dimensions_(vector::zero),
    k_(1),
    s_(vector::zero)
{
    dimensions_ = heatSourceModelCoeffs_.lookup<vector>("dimensions");
    k_ = heatSourceModelCoeffs_.lookup<scalar>("k");

    minimumDepth_ = dimensions_.z();
    depth_ = minimumDepth_;

    D4Sigma_ =
        heatSourceProfiles::superGaussianProfile::D4Sigma
        (
            dimensions_.x(),
            dimensions_.y(),
            k_
        );

    updateSource();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModels::superGaussian::updateSource()
{
    dimensions_.z() = depth_;

    s_ = dimensions_/Foam::pow(2.0, 1.0/k_);

    const scalar fMax =
        Foam::pow
        (
            invIncGammaRatio_P(3.0/k_, 1.0 - tolerance_),
            1.0/k_
        );

    sourceBb_ =
        boundBox
        (
            point(-fMax*s_.x(), -fMax*s_.y(), -fMax*s_.z()),
            point(fMax*s_.x(), fMax*s_.y(), 0)
        );

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            (2.0/3.0)*s_.x()*s_.y()*s_.z()*pi
           *Foam::tgamma(1.0 + 3.0/k_)
        );
}


void Foam::heatSourceModels::superGaussian::update()
{
    updateDepth();
    updateSource();
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
        heatSourceModelCoeffs_.lookup("dimensions") >> dimensions_;
        heatSourceModelCoeffs_.lookup("k") >> k_;

        minimumDepth_ = dimensions_.z();
        depth_ = max(minimumDepth_, isoDepth_);

        D4Sigma_ =
            heatSourceProfiles::superGaussianProfile::D4Sigma
            (
                dimensions_.x(),
                dimensions_.y(),
                k_
            );

        updateSource();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
