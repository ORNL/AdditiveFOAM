/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Oak Ridge National Laboratory
     \\/     M anipulation  |
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

#include "constantAbsorption.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace absorptionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(absorptionModel, constant, dictionary);
}
}

Foam::absorptionModels::constant::constant
(
    const dictionary& dict,
    const fvMesh&
)
:
    eta_("eta", dimless, dict)
{}

// ************************************************************************* //
