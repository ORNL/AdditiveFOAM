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
    Foam::functionObjects::ExaCA

SourceFiles
    ExaCA.C

\*---------------------------------------------------------------------------*/

#ifndef ExaCA_H
#define ExaCA_H

#include "fvMeshFunctionObject.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class ExaCA Declaration
\*---------------------------------------------------------------------------*/

class ExaCA
:
    public fvMeshFunctionObject
{
    // Private Data

        const volScalarField& T_;

        volPointInterpolation vpi_;

        mutable pointScalarField Tp_;

        DynamicList<label> overlapCells;

        DynamicList<List<scalar>> events;

        List<DynamicList<label>> pointsInCell;

        DynamicList<point> positions;
        
        DynamicList<scalarField> weights;
              
        boundBox  box_;

        scalar    isoValue_;

        scalar    dx_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("ExaCA");


    // Constructors

        //- Construct from Time and dictionary
        ExaCA
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ExaCA(const ExaCA&) = delete;


    //- Destructor
    virtual ~ExaCA();


    // Member Functions

        //- Read the ExaCA data
        virtual bool read(const dictionary&);

        //- Initialize lists
        virtual void setOverlapCells();//initialize();

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        virtual void mapPoints();

        virtual void interpolate();
                                
        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the ExaCA data
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ExaCA&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
