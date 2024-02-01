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

#include "absorptionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(absorptionModel, 0);
    defineRunTimeSelectionTable(absorptionModel, dictionary);
}

const Foam::word Foam::absorptionModel::heatSourceDictName
(
    "heatSourceDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::absorptionModel::createIOobject
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
) const
{
    IOobject io
    (
        dict.name(),
        mesh.time().timeName(),
        mesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
    }

    return io;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::absorptionModel::absorptionModel
(
    const word& type,
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(sourceName, dict, mesh)),
    sourceName_(sourceName),
    heatSourceDict_(dict),
    mesh_(mesh),
    
    beamDict_(heatSourceDict_.optionalSubDict(sourceName_)),
    absorptionModelCoeffs_(beamDict_.optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::absorptionModel::read()
{
    if (regIOobject::read())
    {
        absorptionModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
