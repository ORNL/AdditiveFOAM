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
    Foam::refinementModels::refinementVolume

Description
    This refinement control scheme projects along the scan path until the
    global volume of the refined region meets or exceeds a target value
    estimated from a target mesh size.

    The scan-path contribution is counted as the discrete union of newly
    marked cells, so overlapping path vectors and overlapping sources do not
    double count refinement volume.

SourceFiles
    refinementVolume.C

\*---------------------------------------------------------------------------*/

#ifndef refinementVolume_H
#define refinementVolume_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{

/*---------------------------------------------------------------------------*\
                           Class refinementVolume
\*---------------------------------------------------------------------------*/

class refinementVolume
:
    public Foam::refinementModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffs_;


protected:

    // Protected Data

        //- Target cells per processor
        label targetCellsPerProc_;

        //- Target global refinement volume
        scalar targetRefinementVolume_;

        //- Minimum scan-path refinement volume
        scalar minRefinementVolume_;

        //- Time at which to update refinement again
        dimensionedScalar updateTime_;

        //- Ratio of target cells per processor to current cells per processor
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

        //- Update the refinement marker field
        virtual bool update();

        //- Read the heatSourceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace refinementModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
