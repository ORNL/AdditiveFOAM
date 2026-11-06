/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::solidificationData

SourceFiles
    solidificationData.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationData_H
#define solidificationData_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "DynamicList.H"
#include "boundBox.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class solidificationData Declaration
\*---------------------------------------------------------------------------*/

class solidificationData
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Temperature field
        const volScalarField& T_;

        //- Cooling rate field
        volScalarField R_;

        //- Cells that overlap the sampling box
        DynamicList<label> overlapCells;

        //- Solidification events
        DynamicList<List<scalar>> events;

        //- Sampling box
        boundBox box_;

        //- Solidification isotherm
        scalar isoValue_;


    // Private Member Functions

        //- Correct mesh-dependent data
        void correct();


public:

    //- Runtime type information
    TypeName("solidificationData");


    // Constructors

        //- Construct from Time and dictionary
        solidificationData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        solidificationData(const solidificationData&) = delete;


    //- Destructor
    virtual ~solidificationData();


    // Member Functions

        //- Read the solidificationData data
        virtual bool read(const dictionary&);

        //- Initialize lists
        virtual void setOverlapCells();

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute
        virtual bool execute();

        //- Execute at the final time-loop
        virtual bool end();

        //- Write the solidificationData data
        virtual bool write();

        //- Update for mesh motion
        virtual void movePoints(const polyMesh& mesh);

        //- Update for topology change
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update for mesh mapping
        virtual void mapMesh(const polyMeshMap& map);

        //- Update for redistribution
        virtual void distribute(const polyDistributionMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidificationData&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
