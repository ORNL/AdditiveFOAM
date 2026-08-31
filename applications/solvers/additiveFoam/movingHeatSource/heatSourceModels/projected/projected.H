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

Class
    Foam::heatSourceModels::projected

Description
    Separable heat source constructed from a planar beam profile and a
    one-sided axial projection.

SourceFiles
    projected.C

\*---------------------------------------------------------------------------*/

#ifndef projected_H
#define projected_H

#include "heatSourceModel.H"
#include "heatSourceProfile.H"
#include "heatSourceProjection.H"

namespace Foam
{
namespace heatSourceModels
{

class projected
:
    public heatSourceModel
{
    autoPtr<heatSourceProfile> profile_;

    autoPtr<heatSourceProjection> projection_;

    scalar componentTolerance_;

public:

    TypeName("projected");

    projected
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~projected()
    {}

    virtual void update
    (
        const scalar depth,
        const scalar aspectRatio
    );

    virtual scalar weight(const vector& r) const;

    inline virtual const profileMetrics& metrics() const
    {
        return profile_->metrics();
    }
};

} // End namespace heatSourceModels
} // End namespace Foam

#endif

// ************************************************************************* //
