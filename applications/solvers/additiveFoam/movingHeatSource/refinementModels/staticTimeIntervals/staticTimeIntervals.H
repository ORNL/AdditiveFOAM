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
    Foam::refinementModel::staticTimeIntervals

Description
    This refinement control scheme divides the simulation into an even
    number of intervals equal to the simulation end time divided by the
    number of intervals. Refinement for each beam is carried out on each
    interval.

SourceFiles
    staticTimeIntervals.C

\*---------------------------------------------------------------------------*/

#ifndef staticTimeIntervals_H
#define staticTimeIntervals_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModel
{

/*---------------------------------------------------------------------------*\
                            Class staticTimeIntervals
\*---------------------------------------------------------------------------*/

class staticTimeIntervals
:
    public refinementModel
{
    // Private Data
    
        //- Dictionary
        dictionary coeffs_;
        
protected:

    //- Protected data
        
        //- Number of intervals
        scalar intervals_;
        
        //- Time duration of each interval
        scalar intervalLength_;
        
        //- Time for next update
        scalar updateTime_;  

public:

    //- Runtime type information
    TypeName("staticTimeIntervals");


    // Constructors

        //- Construct from components
        staticTimeIntervals
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~staticTimeIntervals()
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
