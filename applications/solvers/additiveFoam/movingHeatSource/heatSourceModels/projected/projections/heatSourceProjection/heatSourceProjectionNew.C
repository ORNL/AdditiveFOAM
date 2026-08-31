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

#include "heatSourceProjection.H"

Foam::autoPtr<Foam::heatSourceProjection> Foam::heatSourceProjection::New
(
    const dictionary& dict
)
{
    const word projectionType(dict.lookup("model"));

    Info<< "Selecting heatSourceProjection " << projectionType << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->find(projectionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << heatSourceProjection::typeName << " type "
            << projectionType << nl << nl
            << "Valid heatSourceProjections are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heatSourceProjection>
    (
        cstrIter()(dict)
    );
}

// ************************************************************************* //
