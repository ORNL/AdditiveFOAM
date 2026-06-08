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

    You should have received a copy of the the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::pathVector

Description
    Normalized representation of a scan path vector with start/end position
    and start/end time. Derived quantities such as displacement, distance,
    duration, and speed are calculated on construction and accessed through
    getter functions.

SourceFiles
    pathVector.C

\*---------------------------------------------------------------------------*/

#ifndef pathVector_H
#define pathVector_H

#include "point.H"
#include "scalar.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pathVector Declaration
\*---------------------------------------------------------------------------*/

class pathVector
{
    // Private Data

        //- Start position of the path vector
        point startPosition_;

        //- End position of the path vector
        point endPosition_;

        //- Start time of the path vector
        scalar startTime_;

        //- End time of the path vector
        scalar endTime_;

        //- Beam power over the path vector
        scalar power_;

        //- Displacement from start position to end position
        vector displacement_;

        //- Magnitude of the displacement
        scalar distance_;

        //- Duration of the path vector
        scalar duration_;

        //- Speed along the path vector
        scalar speed_;


    // Private Member Functions

        //- Set derived quantities
        void setDerived();


public:

    // Constructors

        //- Construct default path vector as a zero-power stationary vector
        pathVector();

        //- Construct from components
        pathVector
        (
            const point& startPosition,
            const point& endPosition,
            const scalar startTime,
            const scalar endTime,
            const scalar power
        );


    // Member Functions

        // Access

            //- Return start position
            inline const point startPosition() const
            {
                return startPosition_;
            }

            //- Return end position
            inline const point endPosition() const
            {
                return endPosition_;
            }

            //- Return start time
            inline scalar startTime() const
            {
                return startTime_;
            }

            //- Return end time
            inline scalar endTime() const
            {
                return endTime_;
            }

            //- Return beam power
            inline scalar power() const
            {
                return power_;
            }

        // Derived

            //- Return displacement
            inline const vector& displacement() const
            {
                return displacement_;
            }

            //- Return distance
            inline scalar distance() const
            {
                return distance_;
            }

            //- Return duration
            inline scalar duration() const
            {
                return duration_;
            }

            //- Return speed
            inline scalar speed() const
            {
                return speed_;
            }

        // Query

            //- Return position at the provided time
            point position(const scalar time) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
