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

Class
    Foam::heatSourceProjection

Description
    Abstract base class for one-dimensional heat source projections.

SourceFiles
    heatSourceProjection.C
    heatSourceProjectionNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatSourceProjection_H
#define heatSourceProjection_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class heatSourceProjection
{
public:

    //- Runtime type information
    TypeName("heatSourceProjection");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        heatSourceProjection,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    // Constructors

        heatSourceProjection()
        {}

        //- Disallow default bitwise copy construction
        heatSourceProjection(const heatSourceProjection&);

    // Selectors

        static autoPtr<heatSourceProjection> New(const dictionary& dict);

    //- Destructor
    virtual ~heatSourceProjection()
    {}

    // Member Functions

        //- Update the projection state
        virtual void update
        (
            const scalar depth,
            const scalar aspectRatio,
            const scalar tolerance
        ) = 0;

        //- Return the projection weight
        virtual scalar weight(const scalar z) const = 0;

        //- Return the projection integral
        virtual scalar integral() const = 0;

        //- Return the lower projection bound
        virtual scalar zMin() const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heatSourceProjection&) = delete;
};

} // End namespace Foam

#endif

// ************************************************************************* //
