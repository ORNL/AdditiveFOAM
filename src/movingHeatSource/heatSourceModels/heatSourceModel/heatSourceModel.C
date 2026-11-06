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

#include "heatSourceModel.H"

namespace Foam
{
    defineTypeNameAndDebug(heatSourceModel, 0);
    defineRunTimeSelectionTable(heatSourceModel, dictionary);
}

Foam::heatSourceModel::heatSourceModel
(
    const dictionary& dict
)
:
    nPoints_(dict.lookupOrDefault<labelVector>("nPoints", vector::one)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 1.0e-3)),
    depth_(dict.lookup<scalar>("depth")),
    sourceDepth_(depth_),
    bounds_(point::zero, point::zero),
    V0_("V0", dimVolume, 0)
{}

// ************************************************************************* //
