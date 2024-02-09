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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::refinementController> Foam::refinementController::New
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
{
    //- Initialize modelType to a non-model word
    word modelType("unselected");
    
    //- Get model type from refinement control subdict
    dictionary refinementControlDict(dict.optionalSubDict("refinementControl"));
    refinementControlDict.lookup("refinementController") >> modelType;

    Info<< "Selecting refinement control model " << modelType << endl;

    //- Look up model type from runtime selection table and throw error
    //  if it doesn't exist
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << refinementController::typeName<< " type "
            << modelType << nl << nl
            << "Valid  refinementControllers are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<refinementController>(cstrIter()(sources, dict, mesh));
}


// ************************************************************************* //
