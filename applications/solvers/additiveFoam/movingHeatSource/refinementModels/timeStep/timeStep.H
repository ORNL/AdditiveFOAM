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
    Foam::refinementModel::timeStep

Description
    This refinement model marks refinementField at every time step, including
    cells which overlap the current beam positions and are above the target
    temperature.

SourceFiles
    timeStep.C

\*---------------------------------------------------------------------------*/

#ifndef timeStep_H
#define timeStep_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{

/*---------------------------------------------------------------------------*\
                            Class timeStep
\*---------------------------------------------------------------------------*/

class timeStep
:
    public refinementModel
{

public:

    //- Runtime type information
    TypeName("timeStep");


    // Constructors

        //- Construct from components
        timeStep
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~timeStep()
    {}


    // Member Functions

        //- Update every time step
        virtual bool update();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace refinementModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
