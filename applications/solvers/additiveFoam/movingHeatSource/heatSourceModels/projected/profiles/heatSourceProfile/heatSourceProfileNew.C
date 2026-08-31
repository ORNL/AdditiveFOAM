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

#include "heatSourceProfile.H"

Foam::autoPtr<Foam::heatSourceProfile> Foam::heatSourceProfile::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const scalar tolerance
)
{
    const word profileType(dict.lookup("model"));

    Info<< "Selecting heatSourceProfile " << profileType << endl;

    const auto cstrIter = dictionaryConstructorTablePtr_->find(profileType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << heatSourceProfile::typeName << " type "
            << profileType << nl << nl
            << "Valid heatSourceProfiles are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<heatSourceProfile>
    (
        cstrIter()(dict, mesh, tolerance)
    );
}

// ************************************************************************* //
