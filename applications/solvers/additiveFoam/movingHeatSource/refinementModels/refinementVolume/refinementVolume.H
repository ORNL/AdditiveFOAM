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
    Foam::refinementModel::refinementVolume

Description
    This refinement control scheme marches along the scan path until the
    volume of the refined region meets or exceeds a target value calculated
    from the target mesh size, and the ratio of the volumes of the refined and
    base-sized cells.

SourceFiles
    refinementVolume.C

\*---------------------------------------------------------------------------*/

#ifndef refinementVolume_H
#define refinementVolume_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModel
{

/*---------------------------------------------------------------------------*\
                           Class refinementVolume
\*---------------------------------------------------------------------------*/

class refinementVolume
:
    public refinementModel
{
    // Private Data
    
        //- Dictionary
        dictionary coeffs_;
        
protected:

    //- Protected data
        
        //- Target cells per processor (default 5k)
        label targetCellsPerProc_;
                
        //- Target volume of refined region
        scalar targetRefinementVolume_;

        //- Minimum refinement volume
        scalar minRefinementVolume_;
        
        //- Time at which to update again
        dimensionedScalar updateTime_;
        
        //- Load imbalance ratio to check AMR performance
        scalar cellLoadBalanceRatio_;

public:

    //- Runtime type information
    TypeName("refinementVolume");


    // Constructors

        //- Construct from components
        refinementVolume
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~refinementVolume()
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
