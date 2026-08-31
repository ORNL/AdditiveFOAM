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

#include "projected.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(projected, 0);
    addToRunTimeSelectionTable
    (
        heatSourceModel,
        projected,
        dictionary
    );
}
}

Foam::heatSourceModels::projected::projected
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(dict),
    profile_
    (
        heatSourceProfile::New
        (
            dict.subDict("profile"),
            mesh,
            tolerance_/(1.0 + Foam::sqrt(1.0 - tolerance_))
        )
    ),
    projection_(heatSourceProjection::New(dict.subDict("projection"))),
    componentTolerance_
    (
        tolerance_/(1.0 + Foam::sqrt(1.0 - tolerance_))
    )
{
    bounds_ = profile_->bounds();

    Info<< "Beam-plane D4Sigma (major minor)="
        << metrics().D4Sigma() << " m, azimuth="
        << metrics().azimuth() << " rad" << endl;
}


void Foam::heatSourceModels::projected::update
(
    const scalar depth,
    const scalar aspectRatio
)
{
    sourceDepth_ = depth;

    projection_->update
    (
        sourceDepth_,
        aspectRatio,
        componentTolerance_
    );

    bounds_.min().z() = projection_->zMin();
    bounds_.max().z() = 0;

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            metrics().integral()*projection_->integral()
        );
}


Foam::scalar Foam::heatSourceModels::projected::weight
(
    const vector& r
) const
{
    return
        projection_->weight(r.z())
       *profile_->weight(r.x(), r.y());
}

// ************************************************************************* //
