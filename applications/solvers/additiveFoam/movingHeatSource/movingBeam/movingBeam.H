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

Class
    Foam::movingBeams

Description
    Class for individual moving beam heat sources

SourceFiles
    movingBeam.C
    calculateQDot.C

\*---------------------------------------------------------------------------*/

#ifndef movingBeam_H
#define movingBeam_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "segment.H"
#include "HashTable.H"
#include "absorptionModel.H"
#include "DynamicList.H"
#include "meshSearch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class movingBeam Declaration
\*---------------------------------------------------------------------------*/

class movingBeam
{
private:

    // Private data
    
        //- Name of beam
        const word& sourceName_;
        
        //- Main moving heat source dictionary
        const dictionary& dict_;
        
        //- Run time information
        const Time& runTime_;
        
        //- Individual beam sub-dict
        const dictionary& beamDict_;

        //- List of beam segments
        DynamicList<segment> path_;
        
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
        
        //- Adjust solution time to hit path segment start times? (d = true)
        bool hitPathIntervals_;

    // Static Data Members

        //- Tolerance for scan path intervals
        static const scalar eps;

public:

    //- Runtime type information
    TypeName("movingBeam");
    
    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        movingBeam,
        dictionary,
        (
            const word& sourceName,
            const dictionary& dict,
            const Time& runTime
        ),
        (sourceName, dict, runTime)
    );

    // Constructors

        // Construct from time and dictionary
        movingBeam
        (
            const word& sourceName,
            const dictionary& dict,
            const Time& runTime
        );
                
        //- Disallow default bitwise copy construction
        movingBeam(const movingBeam&);
        
    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<movingBeam> New
        (
            const word& sourceName,
            const dictionary& dict,
            const Time& runTime
        )
        {
            return autoPtr<movingBeam>(new movingBeam(sourceName, dict, runTime));
        }

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
        
        scalar deltaT() const
        {
            return deltaT_;
        }

        //- Read the path file
        void readPath();
        
        //- Returns true if the simulation time is less than path endTime
        bool activePath();
        
        //- Move the beam to the provided time
        void move(const scalar time);
        
        //- Returns the path index at the provided time
        label findIndex(const scalar time);

        //- Adjust solution time step to hit pathInterval
        void adjustDeltaT(scalar& dt);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

