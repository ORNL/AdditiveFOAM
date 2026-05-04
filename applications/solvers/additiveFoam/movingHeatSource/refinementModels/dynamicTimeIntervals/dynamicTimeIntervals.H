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
    Foam::refinementModel::dynamicTimeIntervals

Description
    This refinement control scheme divides the simulation into sub-intervals
    on which to perform refinement operation. At the start of a new interval,
    the length of the new interval is adjusted to target a specified CPU
    load, in terms of cells per processor. The adjustment is performed using
    the observed CPU load from the previous interval. Relaxation of the
    interval size adjustment is provided to limit oscillations in interval
    and mesh size.

SourceFiles
    dynamicTimeIntervals.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicTimeIntervals_H
#define dynamicTimeIntervals_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModel
{

/*---------------------------------------------------------------------------*\
                           Class dynamicTimeIntervals
\*---------------------------------------------------------------------------*/

class dynamicTimeIntervals
:
    public refinementModel
{
    // Private Data
    
        //- Dictionary
        dictionary coeffs_;
        
protected:

    //- Protected data
        
        //- Target cells per processor (default 10k)
        const int cellsPerProc_;
        
        //- Relaxation factor for interval adjustment (default 0.9)
        const scalar relax_;
        
        //- Characteristic dimension of unrefined base mesh
        scalar unrefinedSize_;
        
        //- Minimum time to maintain at least D4sigma between refinements
        scalar minIntervalTime_;
        
        //- Length of current interval
        scalar intervalLength_;
        
        //- Time at which to update again
        scalar updateTime_;

public:

    //- Runtime type information
    TypeName("dynamicTimeIntervals");


    // Constructors

        //- Construct from components
        dynamicTimeIntervals
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~dynamicTimeIntervals()
    {}


    // Member Functions

        //- Update the refinement marker field and return true at each interval
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
