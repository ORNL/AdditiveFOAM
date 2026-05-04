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
namespace refinementModel
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

Foam::refinementModel::refinementVolume::refinementVolume
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementModel(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    targetCellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 5000)),
    targetRefinementVolume_(0.0),
    minRefinementVolume_(0.0),
    updateTime_(dimTime, 0.0),
    cellLoadBalanceRatio_(1.0)
{
    minRefinementVolume_ = 32.0 * cmptProduct(buffer_);


    Info << "minRefinementVolume_: " << minRefinementVolume_ << endl;

    label nCells0 = mesh_.nCells();
    reduce(nCells0, sumOp<label>());

    // OPTION 1:
    scalar refinementCells_ = minRefinementVolume_ * nCells0 / gSum(mesh_.V());
    
    label minCellsPerCpu_ =
    (
        nCells0
      + refinementCells_ * (Foam::pow(2.0, 3.0 * nLevels_) - 1.0)
    ) / Pstream::nProcs();

    targetCellsPerProc_ = max(targetCellsPerProc_, minCellsPerCpu_);

    targetRefinementVolume_ =
        0.5
      * (gSum(mesh_.V()) / nCells0)
      * max(targetCellsPerProc_*Pstream::nProcs() - nCells0, 0.0)
      / (Foam::pow(2.0, 3.0 * nLevels_) - 1.0);

    Info << "targetRefinementVolume_: " << targetRefinementVolume_ << endl;
    Info << "targetCellsPerProc_: " << targetCellsPerProc_ << endl;
  
    updateTime_ = refinementModel::refineUsingVolume(targetRefinementVolume_);
    
    Info << "Initial AMR update time: " << updateTime_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementModel::refinementVolume::update()
{
    // 0a. Calculate the total number of cells in mesh
    label nCellsTotal = mesh_.nCells();
    reduce(nCellsTotal, sumOp<label>());
    scalar currentCellsPerProc = nCellsTotal / Pstream::nProcs();
    
    // 0b. Calculate cellLoadBalanceRatio from current number of cells and target
    cellLoadBalanceRatio_ = targetCellsPerProc_ / currentCellsPerProc;
       
    if(updateTime_.value() - mesh_.time().value() < small)
    {
        // 1. Reactive refinement based on temperature
        refinementModel::refineUsingTemperature();

        // 2. Scan path completed
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info<< typeName << ": "
                << "Scan path completed. "
                << "Continuing AMR checks for possible mesh coarsening."
                << endl;

            // TODO: Fix hardcoded dilation
            updateTime_ = mesh_.time() + 10*mesh_.time().deltaT();
            return true;
        }

        targetRefinementVolume_ =
            max
            (
                targetRefinementVolume_ * cellLoadBalanceRatio_,
                cmptProduct(buffer_)
            );

        // 3. Predictive refinement using dynamic target volume    
        updateTime_ =
            refinementModel::refineUsingVolume(targetRefinementVolume_, minRefinementVolume_);
            
            
        // TODO: Should we control dilation and shrinking here?
        Info<< typeName << ":" << endl
            << "Target refinement volume: " << targetRefinementVolume_ << endl
            << "Min refinement volume: " << minRefinementVolume_ << endl
            << "Next refinement check is: " << updateTime_ << endl            
            << "Average load imbalance: " << cellLoadBalanceRatio_ << endl;
    }

    return true;
}


bool Foam::refinementModel::refinementVolume::read()
{
    if (refinementModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
