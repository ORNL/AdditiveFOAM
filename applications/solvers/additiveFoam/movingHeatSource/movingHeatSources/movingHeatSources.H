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
    Foam::movingHeatSources

Description
    Aggregate of moving heat sources.

SourceFiles
    movingHeatSources.C

\*---------------------------------------------------------------------------*/

#ifndef movingHeatSources_H
#define movingHeatSources_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "movingHeatSource.H"
#include "refinementModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class movingHeatSources Declaration
\*---------------------------------------------------------------------------*/

class movingHeatSources
{
private:

    // Private data

        //- Pointer to mesh information
        const fvMesh& mesh_;

        //- Main movingHeatSourceDict
        IOdictionary dict_;

        //- Total volumetric power of all beams
        volScalarField qDot_;

        //- Pointer list to each beam instance
        PtrList<movingHeatSource> sources_;

        //- Pointer to refinement model
        autoPtr<refinementModel> refinement_;

        //- Construct sources and refinement model from the dictionary
        void readSources();

public:

    // Constructors

        // Construct from time and dictionary
        movingHeatSources
        (
            const fvMesh& mesh
        );

    //- Destructor
    ~movingHeatSources();

    // Member Functions

        //- Return const access to the total volumetric heating
        const volScalarField& qDot() const
        {
            return qDot_;
        }

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
