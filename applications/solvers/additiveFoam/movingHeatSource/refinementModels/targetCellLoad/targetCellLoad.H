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

    You should have received a copy of the the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::refinementModels::targetCellLoad

Description
    This refinement model estimates a target scan-path refine volume from a
    target cell count per processor. The scan path is projected until the
    global union of marked refinement cells reaches that target volume.

    The target volume is adjusted over time based on the current cell load.

    Example dictionary:

    refinementModel
    {
        refinementModel         targetCellLoad;
        refinementTemperature   1000;

        buffers
        {
            beam                (85e-6 85e-6 100e-6);
        }

        targetCellLoadCoeffs
        {
            targetCellsPerProc          5000;
            nBufferVolumes              4;
            maxSearchIter               10;
            timeTolerance               1e-4;
        }
    }

SourceFiles
    targetCellLoad.C

\*---------------------------------------------------------------------------*/

#ifndef targetCellLoad_H
#define targetCellLoad_H

#include "refinementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{

/*---------------------------------------------------------------------------*\
                           Class targetCellLoad
\*---------------------------------------------------------------------------*/

class targetCellLoad
:
    public Foam::refinementModel
{
    // Private Data

        //- Coefficients dictionary
        dictionary coeffs_;

    // Private Member Functions

        //- Read model coefficients
        void readCoeffs();


protected:

    // Protected Data

        //- Target cells per processor
        label targetCellsPerProc_;

        //- Minimum projected scan-path volume in source-buffer volumes
        scalar nBufferVolumes_;

        //- Maximum number of volume-search bisection iterations
        label maxSearchIter_;

        //- Time tolerance for volume-search bisection
        scalar timeTolerance_;

        //- Fraction of estimated available cell budget used for target volume
        scalar initialTargetVolumeFactor_;

        //- Max growth factor for target refine volume per AMR update
        scalar maxTargetVolumeGrowth_;

        //- Max shrink factor for target refine volume per AMR update
        scalar maxTargetVolumeShrink_;

        //- AMR check interval after scan path completion, in time-step units
        scalar postScanUpdateInterval_;

        //- Target global refine volume
        scalar targetRefineVolume_;

        //- Minimum scan-path refine volume
        scalar minRefineVolume_;

        //- Time at which to update refinement again
        dimensionedScalar updateTime_;

        //- Ratio of target cells per processor to current cells per processor
        scalar targetToCurrentCellRatio_;


public:

    //- Runtime type information
    TypeName("targetCellLoad");


    // Constructors

        //- Construct from components
        targetCellLoad
        (
            const PtrList<heatSourceModel>& sources,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~targetCellLoad()
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
