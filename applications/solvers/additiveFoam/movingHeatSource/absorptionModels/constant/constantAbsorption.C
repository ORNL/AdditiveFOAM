/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace absorptionModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(absorptionModel, constant, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absorptionModels::constant::constant
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionModel(typeName, sourceName, dict, mesh),
    eta_("eta", dimless, absorptionModelCoeffs_)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::absorptionModels::constant::read()
{
    if (absorptionModel::read())
    {
        absorptionModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        absorptionModelCoeffs_.lookup("eta") >> eta_;

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
