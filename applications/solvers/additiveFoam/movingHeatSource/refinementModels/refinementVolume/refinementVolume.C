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

\*---------------------------------------------------------------------------*/

#include "refinementVolume.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModels
{
    defineTypeNameAndDebug(refinementVolume, 0);

    addToRunTimeSelectionTable
    (
        refinementModel,
        refinementVolume,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModels::refinementVolume::refinementVolume
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    Foam::refinementModel(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    targetCellsPerProc_(coeffs_.lookupOrDefault<label>("cellsPerProc", 5000)),
    targetRefinementVolume_(Zero),
    minRefinementVolume_(Zero),
    updateTime_(dimTime, Zero),
    cellLoadBalanceRatio_(1.0)
{
    const scalar minRefinementVolumeFactor =
        refinementDict_.lookupOrDefault<scalar>
        (
            "minRefinementVolumeFactor",
            1.0
        );

    //- Characteristic symmetric spot volume:
    //  (2*buffer.x) * (2*buffer.y) * (2*buffer.z)
    minRefinementVolume_ =
        minRefinementVolumeFactor*8.0*cmptProduct(buffer_);

    Info<< typeName << ": Minimum scan-path refinement volume: "
        << minRefinementVolume_ << endl;

    label nCellsTotal = mesh_.nCells();
    reduce(nCellsTotal, sumOp<label>());

    const scalar meshVolume = gSum(mesh_.V());

    const scalar averageCellVolume = meshVolume/scalar(nCellsTotal);

    const scalar refinementCellMultiplier =
        Foam::pow(2.0, 3.0*nLevels_) - 1.0;

    if (refinementCellMultiplier <= small)
    {
        FatalErrorInFunction
            << typeName << " requires nLevels > 0"
            << nl
            << exit(FatalError);
    }

    const scalar minimumRefinementCells =
        minRefinementVolume_/averageCellVolume*refinementCellMultiplier;

    const scalar minimumCellsPerProcScalar =
        (scalar(nCellsTotal) + minimumRefinementCells)/Pstream::nProcs();

    label minimumCellsPerProc = label(minimumCellsPerProcScalar);

    if (scalar(minimumCellsPerProc) < minimumCellsPerProcScalar)
    {
        ++minimumCellsPerProc;
    }

    targetCellsPerProc_ =
        max(targetCellsPerProc_, minimumCellsPerProc);

    const scalar targetTotalCells =
        scalar(targetCellsPerProc_)*Pstream::nProcs();

    const scalar availableRefinedCells =
        max(targetTotalCells - scalar(nCellsTotal), scalar(0));

    const scalar targetVolumeSafetyFactor =
        coeffs_.lookupOrDefault<scalar>
        (
            "targetVolumeSafetyFactor",
            0.5
        );

    targetRefinementVolume_ =
        max
        (
            targetVolumeSafetyFactor
           *averageCellVolume
           *availableRefinedCells
           /refinementCellMultiplier,
            minRefinementVolume_
        );

    Info<< typeName << ": Target refinement volume: "
        << targetRefinementVolume_ << endl;

    Info<< typeName << ": Target cells per processor: "
        << targetCellsPerProc_ << endl;

    updateTime_ =
        Foam::refinementModel::refineUsingVolume
        (
            dimensionedScalar
            (
                "targetRefinementVolume",
                dimVolume,
                targetRefinementVolume_
            )
        );

    Info<< typeName << ": Initial AMR update time: "
        << updateTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementModels::refinementVolume::update()
{
    label nCellsTotal = mesh_.nCells();
    reduce(nCellsTotal, sumOp<label>());

    const scalar currentCellsPerProc =
        scalar(nCellsTotal)/Pstream::nProcs();

    cellLoadBalanceRatio_ =
        scalar(targetCellsPerProc_)/max(currentCellsPerProc, small);

    if ((updateTime_.value() - mesh_.time().value()) < small)
    {
        //- Reactive refinement based on temperature
        Foam::refinementModel::refineUsingTemperature();

        //- Scan path completed
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info<< typeName << ": "
                << "Scan path completed. "
                << "Continuing AMR checks for possible mesh coarsening."
                << endl;

            const scalar postScanUpdateIntervalFactor =
                coeffs_.lookupOrDefault<scalar>
                (
                    "postScanUpdateIntervalFactor",
                    10.0
                );

            updateTime_ =
                mesh_.time()
              + postScanUpdateIntervalFactor*mesh_.time().deltaT();

            return true;
        }

        const scalar maxTargetVolumeGrowth =
            coeffs_.lookupOrDefault<scalar>
            (
                "maxTargetVolumeGrowth",
                1.25
            );

        const scalar maxTargetVolumeShrink =
            coeffs_.lookupOrDefault<scalar>
            (
                "maxTargetVolumeShrink",
                0.8
            );

        const scalar targetVolumeRatio =
            min
            (
                max(cellLoadBalanceRatio_, maxTargetVolumeShrink),
                maxTargetVolumeGrowth
            );

        targetRefinementVolume_ =
            max
            (
                targetRefinementVolume_*targetVolumeRatio,
                minRefinementVolume_
            );

        updateTime_ =
            Foam::refinementModel::refineUsingVolume
            (
                dimensionedScalar
                (
                    "targetRefinementVolume",
                    dimVolume,
                    targetRefinementVolume_
                )
            );

        Info<< typeName << ":" << endl
            << "    Target refinement volume: "
            << targetRefinementVolume_ << endl
            << "    Minimum scan-path refinement volume: "
            << minRefinementVolume_ << endl
            << "    Next refinement check: "
            << updateTime_ << endl
            << "    Target/current cell ratio: "
            << cellLoadBalanceRatio_ << endl;
    }

    return true;
}


bool Foam::refinementModels::refinementVolume::read()
{
    return Foam::refinementModel::read();
}


// ************************************************************************* //
