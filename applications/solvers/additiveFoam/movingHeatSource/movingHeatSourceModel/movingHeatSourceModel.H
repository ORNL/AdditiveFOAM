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
    Foam::movingHeatSourceModel

Description
    Class for multiple moving beam heat sources

SourceFiles
    movingHeatSourceModel.C
    calculateQDot.C

\*---------------------------------------------------------------------------*/

#ifndef movingHeatSourceModel_H
#define movingHeatSourceModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "heatSourceModel.H"
#include "movingBeam.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class movingHeatSourceModel Declaration
\*---------------------------------------------------------------------------*/

class movingHeatSourceModel
{
private:

    // Private data

        //- Pointer to mesh information
        const fvMesh& mesh_;
        
        //- Main movingHeatSourceDict
        IOdictionary dict_;
        
        //- List of individual heat sources
        wordList sourceNames_;
        
        //- Subdictionary for refinement controls
        dictionary refineDict_;
        
        //- Switch for turning on refinement field
        bool refine_;
        
        //- Switch for controlling refinement in melt pool behind beam
        bool resolveTail_;
        
        //- Number of time steps between refinements
        label refinementInterval_;
        
        //- Max refinement level specified in dynamicMeshDict. If > 1,
        //  mesh will be updated in successive time steps to ensure cells
        //  newly marked for refinement are completely refined. This overcomes
        //  the limitation of only being able to refine one level per mesh update
        label refinementLevel_;
        
        //- Counter variable for next time index to refine
        label refinementIndex_;
        
        //- How many time steps ahead of the beam to integrate. By setting
        //  fwdSteps = n * refineInterval, cells ahead of the beam will go through
        //  n refinements before the beam reaches their location. This overcomes
        //  the limitation of the refinement function, which can only refine
        //  one level per function call, and ensures the cells are refined
        //  to the desired level before the beam reaches their location. This is an
        //  alternative strategy to the refinementLevel_ control, which
        //  refines n successive times in a row once the refinement interval is reached.
        label fwdSteps_;
        
        //- How many time steps behind beam to integrate. This is generally mutually 
        //  exclusive with the resolveTail_ control, and may be more stable than
        //  refining based on the phase fraction
        label revSteps_;
        
        //- Total volumetric power of all beams
        volScalarField qDot_;
        
        //- Refinement marker field
        volScalarField refinementField_;
                
        //- Pointer list to each beam instance
        PtrList<heatSourceModel> sources_;

public:

    // Constructors

        // Construct from time and dictionary
        movingHeatSourceModel
        (
            const fvMesh& mesh
        );

    //- Destructor
    ~movingHeatSourceModel();

    // Member Functions
    
        //- Return const access to the total volumetric heating
        const volScalarField& qDot() const
        {
            return qDot_;
        }
        
        //- Return const access to the refinement marker field
        const volScalarField& refinementField() const
        {
            return refinementField_;
        }
        
        //- Return true if it is time to refine
        bool refineTime();
        
        //- Adjust deltaT using the current state of each beam
        void adjustDeltaT(scalar& deltaT);
        
        //- Update total qDot field
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

