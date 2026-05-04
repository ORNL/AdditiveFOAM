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

#include "dynamicTimeIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementModel
{
    defineTypeNameAndDebug(dynamicTimeIntervals, 0);
    addToRunTimeSelectionTable
    (
        refinementModel,
        dynamicTimeIntervals,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModel::dynamicTimeIntervals::dynamicTimeIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementModel(typeName, sources, dict, mesh),
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    cellsPerProc_(coeffs_.lookupOrDefault<int>("cellsPerProc", 10000)),
    relax_(coeffs_.lookupOrDefault<scalar>("relax", 0.9)),
    unrefinedSize_(coeffs_.lookupOrDefault<scalar>("unrefinedSize", -1.0)),
    minIntervalTime_(0.0),
    intervalLength_(0.0),
    updateTime_(0.0)
{
    //- Search beams to find maximum number of intervals, corresponding to
    //  the beam with the highest ratio of scan path length to beam width
    scalar maxIntervals = 0.0;

    forAll(sources_, i)
    {
        const scalar bbLen = sources_[i].beam().totalLength();
        const label bbSpots = sources_[i].beam().nSpots();
        
        treeBoundBox beamBb
        (
            min(-1.5 * sources_[i].dimensions(), -buffer_),
            max(1.5 * sources_[i].dimensions(), buffer_)
        );
        
        point bbMin = beamBb.min();
        point bbMax = beamBb.max();
        
        const scalar bbMaxDim = max(bbMax[0] - bbMin[0], bbMax[1] - bbMin[1]);
        
        maxIntervals = max(maxIntervals, bbLen / bbMaxDim + bbSpots);
    }
    
    //- Set minimum interval time from the maximum number of intervals
    minIntervalTime_ = endTime_ / maxIntervals;
    
    //- Target mesh size
    scalar targetCells = cellsPerProc_ * Pstream::nProcs();
    
    //- Find initial mesh size
    scalar nCells0 = mesh_.nCells();
    reduce(nCells0, sumOp<scalar>());
    
    //- Provide warning if target mesh size is smaller than initial mesh size
    if (nCells0 > targetCells)
    {
        Info << "dynamicTimeIntervals: WARNING - initial mesh size larger than "
                "target mesh size." << endl;
    }
    //- Otherwise, estimate refined volume required to hit target mesh size
    else
    {
        //- Estimate unrefined mesh size if not provided
        if (unrefinedSize_ < 0.0)
        {
            unrefinedSize_ = gSum(mesh_.V()) / nCells0;
            
            Info << "dynamicTimeIntervals: estimated unrefined mesh size "
                 << "to be " << unrefinedSize_ << " m." << endl;
        }
        
        //- Estimate volume of refined scan path to hit target mesh size
        const scalar refVol
            = Foam::pow(unrefinedSize_, 3.0) * (targetCells - nCells0)
              / (Foam::pow(2.0, 3.0 * nLevels_) - 1.0);
              
        //- Refine first interval using volume estimate, and set next update
        //  time using the refinementModel::refineUsingVolume function
        updateTime_ =
            refinementModel::refineUsingVolume
            (
                refVol,
                minIntervalTime_
            ).value();
        
        intervalLength_
            = max(updateTime_ - mesh_.time().value(), minIntervalTime_);
            
        const scalar intervals = endTime_ / intervalLength_;
        
        Info << "dynamicTimeIntervals: estimated length of first interval to "
             << "be " << intervalLength_ << " s from refined scan path volume,"
             << " corresponding to " << intervals << " intervals." << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementModel::dynamicTimeIntervals::update()
{
    //- Update if mesh time equals update time
    //  OR if time index is equal to the max refinement level.
    //  This second condition adjusts the mesh after the guess at the first
    //  refinement interval size to prevent an overly long first interval.
    if ((updateTime_ - mesh_.time().value() < small)
      ||
        (mesh_.time().timeIndex() == nLevels_ + 1))
    {
        //- Refine in regions above specified temperature
        refinementModel::refineUsingTemperature();

        //- Don't perform additional refinements if scan path is completed
        if ((endTime_ - mesh_.time().value()) < small)
        {
            Info << "dynamicTimeIntervals: Scan path completed. Continuing AMR"
                 << " checks for possible mesh coarsening" << endl;
                 
            updateTime_ = mesh_.time().value() + intervalLength_;
                 
            return true;
        }
        
        //- Calculate current CPU load (cells per processor)
        label nCells = mesh_.nCells();
        reduce(nCells, sumOp<label>());
        scalar currCellsPerProc = nCells / Pstream::nProcs();
        
        Info << "dynamicTimeIntervals: Current CPU load is "
             << currCellsPerProc
             << " cells per processor, with an interval of length "
             << intervalLength_ << " s." << endl;
        
        //- Rescale interval length
        intervalLength_ =
            relax_ * cellsPerProc_ / currCellsPerProc * intervalLength_
          + (1.0 - relax_) * intervalLength_;
              
        intervalLength_ = max(intervalLength_, minIntervalTime_);
        
        //- Update next refinement time
        updateTime_ = mesh_.time().value() + intervalLength_;
        
        Info << "dynamicTimeIntervals: rescaled interval to "
             << intervalLength_ << " s. Next update will occur at " 
             << updateTime_ << "s. Updating AMR marker field." << endl;
        
        //- Update marker field using calculated update time
        refinementModel::refineUsingTime(updateTime_);
    }

    return true;
}


bool Foam::refinementModel::dynamicTimeIntervals::read()
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
