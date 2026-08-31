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

Class
    Foam::heatSourceModels::modifiedSuperGaussian

Description
    One-sided modified super-Gaussian heat source with finite depth.

SourceFiles
    modifiedSuperGaussian.C

\*---------------------------------------------------------------------------*/

#ifndef modifiedSuperGaussian_H
#define modifiedSuperGaussian_H

#include "heatSourceModel.H"

namespace Foam
{
namespace heatSourceModels
{

class modifiedSuperGaussian
:
    public heatSourceModel
{
    vector2D radius_;
    scalar k_;
    scalar m_;
    scalar coefficient_;
    scalar cosTheta_;
    scalar sinTheta_;
    profileMetrics metrics_;

public:

    TypeName("modifiedSuperGaussian");

    modifiedSuperGaussian
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~modifiedSuperGaussian()
    {}

    virtual void update
    (
        const scalar depth,
        const scalar aspectRatio
    );

    virtual scalar weight(const vector& r) const;

    inline virtual const profileMetrics& metrics() const
    {
        return metrics_;
    }
};

} // End namespace heatSourceModels
} // End namespace Foam

#endif

// ************************************************************************* //
