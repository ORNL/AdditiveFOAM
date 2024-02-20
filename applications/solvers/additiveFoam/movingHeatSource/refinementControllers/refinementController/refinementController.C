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

#include "refinementController.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinementController, 0);
    defineRunTimeSelectionTable(refinementController, dictionary);
}

const Foam::word Foam::refinementController::refinementControllerDictName
(
    "refinementControllerDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::refinementController::createIOobject
(
    const dictionary& dict,
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        dict.name(),
        mesh.time().constant(),
        mesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementController::refinementController
(
    const word& type,
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(dict, mesh)),

    sources_(sources),
    mesh_(mesh),
    heatSourceDict_(dict),
    refinementDict_(heatSourceDict_.optionalSubDict("refinementControl")),
    refine_
    (
        (type != "none")
      ? refinementDict_.lookup<bool>("refine")
      : false
    ),
    nLevels_
    (
        (type != "none")
      ? refinementDict_.lookup<label>("nLevels")
      : 0
    ),
    refinementTemperature_
    (
        (type != "none")
      ? refinementDict_.lookupOrDefault<scalar>("refinementTemperature", GREAT)
      : GREAT
    ),
    lastRefinementIndex_(0),
    refinementField_
    (
        IOobject
        (
            "refinementField",
            mesh_.time().timeName(),
            mesh_,
            refine_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            refine_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    )
{
    Info << "refinement temperature: " << refinementTemperature_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementController::update()
{
    if (mesh_.time().timeIndex() == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void Foam::refinementController::setRefinementField()
{
    // TODO: Add gradient based criteria
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    forAll(mesh_.cells(), celli)
    {
        if (T[celli] >= refinementTemperature_)
        {
            refinementField_[celli] = 1;
        }
        else
        {
            refinementField_[celli] = 0;
        }
    }

    refinementField_.correctBoundaryConditions();
}

bool Foam::refinementController::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
