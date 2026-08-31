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
    Foam::heatSourceProfiles::superGaussian

Description
    Planar elliptical super-Gaussian intensity profile. radius contains the
    local x/y radii and need not be ordered. definition is e2 or secondMoment.

SourceFiles
    superGaussianProfile.C

\*---------------------------------------------------------------------------*/

#ifndef superGaussianProfile_H
#define superGaussianProfile_H

#include "heatSourceProfile.H"

namespace Foam
{
namespace heatSourceProfiles
{

class superGaussian
:
    public heatSourceProfile
{
    vector2D radius_;
    scalar k_;
    scalar coefficient_;
    scalar cosTheta_;
    scalar sinTheta_;
    profileMetrics metrics_;
    boundBox bounds_;

public:

    TypeName("superGaussian");

    superGaussian
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const scalar tolerance
    );

    static scalar coefficient
    (
        const dictionary& dict,
        const scalar k
    );

    virtual ~superGaussian()
    {}

    virtual scalar weight(const scalar x, const scalar y) const;

    inline virtual const profileMetrics& metrics() const
    {
        return metrics_;
    }

    inline virtual const boundBox& bounds() const
    {
        return bounds_;
    }
};

} // End namespace heatSourceProfiles
} // End namespace Foam

#endif

// ************************************************************************* //
