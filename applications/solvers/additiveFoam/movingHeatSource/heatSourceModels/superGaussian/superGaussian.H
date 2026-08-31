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
    Foam::heatSourceModels::superGaussian

Description
    One-sided volumetric super-Gaussian heat source. definition applies to
    radius; depth always specifies the axial depth.

SourceFiles
    superGaussian.C

\*---------------------------------------------------------------------------*/

#ifndef superGaussian_H
#define superGaussian_H

#include "heatSourceModel.H"

namespace Foam
{
namespace heatSourceModels
{

class superGaussian
:
    public heatSourceModel
{
    vector2D radius_;

    scalar k_;

    scalar coefficient_;

    scalar cosTheta_;

    scalar sinTheta_;

    profileMetrics metrics_;

public:

    TypeName("superGaussian");

    superGaussian
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~superGaussian()
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
