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

\*---------------------------------------------------------------------------*/

#include "targetCellLoad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{
    defineTypeNameAndDebug(targetCellLoad, 0);

    addToRunTimeSelectionTable
    (
        refinementModel,
        targetCellLoad,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementModels::targetCellLoad::readCoeffs()
{
    targetCellsPerProc_ =
        coeffs_.lookupOrDefault<label>("targetCellsPerProc", 5000);

    nBufferVolumes_ =
        coeffs_.lookupOrDefault<scalar>
        (
            "nBufferVolumes",
            1.0
        );

    maxSearchIter_ =
        coeffs_.lookupOrDefault<label>("maxSearchIter", 30);

    timeTolerance_ =
        coeffs_.lookupOrDefault<scalar>("timeTolerance", small);

    initialTargetVolumeFactor_ =
        coeffs_.lookupOrDefault<scalar>
        (
            "initialTargetVolumeFactor",
            0.5
        );

    maxTargetVolumeGrowth_ =
        coeffs_.lookupOrDefault<scalar>("maxTargetVolumeGrowth", 1.2);

    maxTargetVolumeShrink_ =
        coeffs_.lookupOrDefault<scalar>("maxTargetVolumeShrink", 0.8);

    postScanUpdateInterval_ =
        coeffs_.lookupOrDefault<scalar>
        (
            "postScanUpdateInterval",
            10.0
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModels::targetCellLoad::targetCellLoad
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    Foam::refinementModel(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    targetCellsPerProc_
    (
        coeffs_.lookupOrDefault<label>("targetCellsPerProc", 5000)
    ),
    nBufferVolumes_
    (
        coeffs_.lookupOrDefault<scalar>
        (
            "nBufferVolumes",
            1.0
        )
    ),
    maxSearchIter_
    (
        coeffs_.lookupOrDefault<label>("maxSearchIter", 30)
    ),
    timeTolerance_
    (
        coeffs_.lookupOrDefault<scalar>("timeTolerance", small)
    ),
    initialTargetVolumeFactor_
    (
        coeffs_.lookupOrDefault<scalar>
        (
            "initialTargetVolumeFactor",
            0.5
        )
    ),
    maxTargetVolumeGrowth_
    (
        coeffs_.lookupOrDefault<scalar>("maxTargetVolumeGrowth", 1.2)
    ),
    maxTargetVolumeShrink_
    (
        coeffs_.lookupOrDefault<scalar>("maxTargetVolumeShrink", 0.8)
    ),
    postScanUpdateInterval_
    (
        coeffs_.lookupOrDefault<scalar>
        (
            "postScanUpdateInterval",
            10.0
        )
    ),
    targetRefineVolume_(Zero),
    minRefineVolume_(Zero),
    updateTime_(dimTime, Zero),
    targetToCurrentCellRatio_(1.0)
{
    minRefineVolume_ =
        Foam::refinementModel::minimumScanPathRefineVolume
        (
            nBufferVolumes_
        );

    Info<< typeName << ": Minimum scan-path refine volume: "
        << minRefineVolume_ << endl;

    label nCellsTotal = mesh_.nCells();
    reduce(nCellsTotal, sumOp<label>());

    const scalar meshVolume = gSum(mesh_.V());

    const scalar averageCellVolume = meshVolume/nCellsTotal;

    const scalar refineCellMultiplier =
        max(Foam::pow(2.0, 3.0*maxRefinementLevel_) - 1.0, 1.0);

    const scalar minRefineCells =
        minRefineVolume_/averageCellVolume*refineCellMultiplier;

    const label minCellsPerProc =
        (nCellsTotal + minRefineCells)/Pstream::nProcs();

    targetCellsPerProc_ =
        max(targetCellsPerProc_, minCellsPerProc);

    const scalar targetTotalCells =
        targetCellsPerProc_*Pstream::nProcs();

    const scalar availableRefineCells =
        max(targetTotalCells - nCellsTotal, 0.0);

    targetRefineVolume_ =
        max
        (
            initialTargetVolumeFactor_
           *averageCellVolume
           *availableRefineCells
           /refineCellMultiplier,
            minRefineVolume_
        );

    Info<< typeName << ": Target refine volume: "
        << targetRefineVolume_ << endl;

    Info<< typeName << ": Target cells per processor: "
        << targetCellsPerProc_ << endl;

    updateTime_ =
        Foam::refinementModel::markScanPathVolume
        (
            dimensionedScalar
            (
                "targetRefineVolume",
                dimVolume,
                targetRefineVolume_
            ),
            minRefineVolume_,
            maxSearchIter_,
            timeTolerance_
        );

    Info<< typeName << ": Initial AMR update time: "
        << updateTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementModels::targetCellLoad::update()
{
    label nCellsTotal = mesh_.nCells();
    reduce(nCellsTotal, sumOp<label>());

    const scalar currentCellsPerProc =
        scalar(nCellsTotal)/Pstream::nProcs();

    targetToCurrentCellRatio_ =
        targetCellsPerProc_/currentCellsPerProc;

    if ((updateTime_.value() - mesh_.time().value()) < small)
    {
        Foam::refinementModel::markTemperature();

        if ((scanEndTime_ - mesh_.time().value()) < small)
        {
            Info<< typeName << ": "
                << "Scan path completed. "
                << "Continuing AMR checks for possible mesh coarsening."
                << endl;

            updateTime_ =
                mesh_.time()
              + postScanUpdateInterval_*mesh_.time().deltaT();

            return true;
        }

        const scalar targetVolumeRatio =
            min
            (
                max(targetToCurrentCellRatio_, maxTargetVolumeShrink_),
                maxTargetVolumeGrowth_
            );

        targetRefineVolume_ =
            max
            (
                targetRefineVolume_*targetVolumeRatio,
                minRefineVolume_
            );

        updateTime_ =
            Foam::refinementModel::markScanPathVolume
            (
                dimensionedScalar
                (
                    "targetRefineVolume",
                    dimVolume,
                    targetRefineVolume_
                ),
                minRefineVolume_,
                maxSearchIter_,
                timeTolerance_
            );

        Info<< typeName << ":" << endl
            << "    Target refine volume: "
            << targetRefineVolume_ << endl
            << "    Minimum scan-path refine volume: "
            << minRefineVolume_ << endl
            << "    Next refinement check: "
            << updateTime_ << endl
            << "    Target/current cell ratio: "
            << targetToCurrentCellRatio_ << endl;
    }

    return true;
}


bool Foam::refinementModels::targetCellLoad::read()
{
    if (Foam::refinementModel::read())
    {
        coeffs_ = refinementDict_.optionalSubDict(typeName + "Coeffs");

        readCoeffs();

        minRefineVolume_ =
            Foam::refinementModel::minimumScanPathRefineVolume
            (
                nBufferVolumes_
            );

        return true;
    }

    return false;
}


// ************************************************************************* //
