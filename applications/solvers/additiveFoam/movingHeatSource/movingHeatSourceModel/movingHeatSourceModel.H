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

#include "fvCFD.H"
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
        
        //- Total volumetric power of all beams
        volScalarField qDot_;
                
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
        
        //- Adjust deltaT for each beam
        void adjustDeltaT(scalar& deltaT);
        
        //- Return boolean for whether a heat source has transient dimensions
        bool transientExists();

        //- Update total qDot field
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

