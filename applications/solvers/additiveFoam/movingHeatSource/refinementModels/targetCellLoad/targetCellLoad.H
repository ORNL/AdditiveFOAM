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

Class
    Foam::refinementModels::targetCellLoad

Description
    This refinement model estimates a target scan-path refine volume from a
    target cell count per processor. The scan path is projected until the
    global union of marked refinement cells reaches that target volume.

    The target volume is adjusted over time based on the current cell load.

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


protected:

    // Protected Data

        //- Target cells per processor
        label targetCellsPerProc_;

        //- Fraction of estimated available cell budget used for target volume
        scalar initialTargetVolumeFactor_;

        //- Max growth factor for target refine volume per AMR update
        scalar maxTargetVolumeGrowth_;

        //- Max shrink factor for target refine volume per AMR update
        scalar maxTargetVolumeShrink_;

        //- AMR check interval after scan path completion, in time-step units
        scalar postScanUpdateIntervalFactor_;

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
