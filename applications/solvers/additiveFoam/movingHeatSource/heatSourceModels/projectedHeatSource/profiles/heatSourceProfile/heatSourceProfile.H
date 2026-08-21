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
    Foam::heatSourceProfile

Description
    Abstract base class for planar heat source profiles.

SourceFiles
    heatSourceProfile.C
    heatSourceProfileNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatSourceProfile_H
#define heatSourceProfile_H

#include "autoPtr.H"
#include "boundBox.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class fvMesh;

class heatSourceProfile
{
public:

    //- Runtime type information
    TypeName("heatSourceProfile");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        heatSourceProfile,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const scalar tolerance
        ),
        (dict, mesh, tolerance)
    );

    // Constructors

        heatSourceProfile()
        {}

        //- Disallow default bitwise copy construction
        heatSourceProfile(const heatSourceProfile&);

    // Selectors

        static autoPtr<heatSourceProfile> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const scalar tolerance
        );

    //- Destructor
    virtual ~heatSourceProfile()
    {}

    // Member Functions

        //- Return the profile weight
        virtual scalar weight(const scalar x, const scalar y) const = 0;

        //- Return the planar profile integral
        virtual scalar integral() const = 0;

        //- Return the area-equivalent second-moment diameter
        virtual scalar D4Sigma() const = 0;

        //- Return the planar profile bounds
        virtual const boundBox& profileBb() const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const heatSourceProfile&) = delete;
};

} // End namespace Foam

#endif

// ************************************************************************* //
