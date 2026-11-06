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
    Foam::movingBeam

Description
    Class for individual moving beam heat sources

SourceFiles
    movingBeam.C

\*---------------------------------------------------------------------------*/

#ifndef movingBeam_H
#define movingBeam_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "DynamicList.H"
#include "pathVector.H"
#include "Time.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class movingBeam Declaration
\*---------------------------------------------------------------------------*/

class movingBeam
{
private:

    // Private data

        //- Run time information
        const Time& runTime_;

        //- List of beam pathVectors
        DynamicList<pathVector> path_;

        //- Index of path
        label index_;

        //- Position of beam
        vector position_;

        //- Power of beam
        scalar power_;

        //- End time of path
        scalar endTime_;

        //- Time step used for heat source integration
        scalar deltaT_;

        //- Adjust solution time to hit path vector start times?
        bool hitPathIntervals_;

    // Static Data Members

        //- Tolerance for scan path intervals
        static const scalar eps;

public:

    //- Runtime type information
    TypeName("movingBeam");

    // Constructors

        // Construct from time and dictionary
        movingBeam
        (
            const dictionary& sourceDict,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        movingBeam(const movingBeam&);

    //- Destructor
    virtual ~movingBeam()
    {}

    // Member Functions

        //- Return current position of the moving beam
        inline vector position() const
        {
            return position_;
        }

        //- Return current power of the moving beam
        inline scalar power() const
        {
            return power_;
        }

        //- Return time step for path
        inline scalar deltaT() const
        {
            return deltaT_;
        }

        //- Return path end time
        inline scalar endTime() const
        {
            return endTime_;
        }

        //- Return access to the scan path vectors
        inline const DynamicList<pathVector>& path() const
        {
            return path_;
        }

         //- Returns the path index at the provided time
        label getIndex(const scalar time) const;

        //- Returns true if the simulation time is less than path endTime
        bool activePath() const;

        //- Adjust solution time step to hit pathInterval
        void adjustDeltaT(scalar& dt) const;

        //- Move the beam to the provided time
        void move(const scalar time);

        //- Read the path file
        void readPath(const fileName& pathName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
