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

#include "projectedHeatSource.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(projectedHeatSource, 0);
    addToRunTimeSelectionTable
    (
        heatSourceModel,
        projectedHeatSource,
        dictionary
    );
}
}

Foam::heatSourceModels::projectedHeatSource::projectedHeatSource
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    profile_(nullptr),
    projection_(nullptr),
    componentTolerance_
    (
        // (1 - componentTolerance)^2 = 1 - tolerance
        tolerance_/(1.0 + Foam::sqrt(1.0 - tolerance_))
    )
{
    profile_ = heatSourceProfile::New
    (
        heatSourceModelCoeffs_,
        mesh_,
        componentTolerance_
    );

    projection_ = heatSourceProjection::New(heatSourceModelCoeffs_);

    minimumDepth_ = heatSourceModelCoeffs_.lookup<scalar>("minimumDepth");

    depth_ = minimumDepth_;
    D4Sigma_ = profile_->D4Sigma();
    sourceBb_ = profile_->profileBb();

    updateSource();

    Info<< "D4Sigma: " << D4Sigma_ << " m" << endl;
}

void Foam::heatSourceModels::projectedHeatSource::updateSource()
{
    projection_->update
    (
        depth_,
        depth_/(0.5*D4Sigma_),
        componentTolerance_
    );

    sourceBb_.min().z() = projection_->zMin();
    sourceBb_.max().z() = 0;

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            profile_->integral()*projection_->integral()
        );
}

void Foam::heatSourceModels::projectedHeatSource::update()
{
    updateDepth();
    updateSource();
}

inline Foam::scalar
Foam::heatSourceModels::projectedHeatSource::weight(const vector& r) const
{
    const scalar p = projection_->weight(r.z());

    return p > 0 ? p*profile_->weight(r.x(), r.y()) : 0;
}

bool Foam::heatSourceModels::projectedHeatSource::read()
{
    if (heatSourceModel::read())
    {
        // (1 - componentTolerance)^2 = 1 - tolerance
        componentTolerance_ =
            tolerance_/(1.0 + Foam::sqrt(1.0 - tolerance_));

        profile_ = heatSourceProfile::New
        (
            heatSourceModelCoeffs_,
            mesh_,
            componentTolerance_
        );

        projection_ = heatSourceProjection::New(heatSourceModelCoeffs_);

        minimumDepth_ = heatSourceModelCoeffs_.lookup<scalar>("minimumDepth");

        depth_ = max(minimumDepth_, isoDepth_);
        D4Sigma_ = profile_->D4Sigma();
        sourceBb_ = profile_->profileBb();

        updateSource();

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
