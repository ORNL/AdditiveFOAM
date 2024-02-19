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
    Foam::refinementControllers::uniformIntervals

Description
    This refinement control scheme divides the simulation into an even
    number of intervals equal to the simulation end time divided by the
    number of intervals. Refinement for each beam is carried out on each
    interval.

SourceFiles
    uniformIntervals.C

\*---------------------------------------------------------------------------*/

#ifndef uniformIntervals_H
#define uniformIntervals_H

#include "refinementController.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{

/*---------------------------------------------------------------------------*\
                             Class uniformIntervals
\*---------------------------------------------------------------------------*/

class uniformIntervals
:
    public refinementController
{
    // Private Data
    
        //- Dictionary
        dictionary coeffs_;
        
        //- Number of intervals
        label intervals_;
        
        //- User-defined bounding box for refinement around moving beam
        //- Defaults to an inverted bounding box that is updated for each beam
        boundBox boundingBox_;
        
        //- Time duration of each interval
        scalar intervalTime_;
        
        //- Time for next update
        scalar updateTime_;

public:

    //- Runtime type information
    TypeName("uniformIntervals");


    // Constructors

        //- Construct from components
        uniformIntervals
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~uniformIntervals()
    {}


    // Member Functions

        //- Update the refinement marker field and return true at each interval
        virtual bool update();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace refinementControllers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
