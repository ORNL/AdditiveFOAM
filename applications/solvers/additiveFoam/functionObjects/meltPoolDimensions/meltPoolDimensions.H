/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

Class
    Foam::functionObjects::meltPoolDimensions

SourceFiles
    meltPoolDimensions.C

\*---------------------------------------------------------------------------*/

#ifndef meltPoolDimensions_H
#define meltPoolDimensions_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "OFstream.H"
#include "PtrList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class meltPoolDimensions Declaration
\*---------------------------------------------------------------------------*/

class meltPoolDimensions
:
    public fvMeshFunctionObject
{
    // Private Data

        const volScalarField& T_;
        scalarList isoValues_;
        scalar scanPathAngle_;

        PtrList<OFstream> logFilePtrs_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("meltPoolDimensions");


    // Constructors

        //- Construct from Time and dictionary
        meltPoolDimensions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        meltPoolDimensions(const meltPoolDimensions&) = delete;


    //- Destructor
    virtual ~meltPoolDimensions();


    // Member Functions

        //- Read the meltPoolDimensions data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();
                                
        //- Execute at the final time-loop
        virtual bool end();

        //- Write the meltPoolDimensions data
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meltPoolDimensions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
