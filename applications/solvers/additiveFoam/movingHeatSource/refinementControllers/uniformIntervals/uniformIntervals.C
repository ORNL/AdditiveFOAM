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

#include "uniformIntervals.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace refinementControllers
{
    defineTypeNameAndDebug(uniformIntervals, 0);
    addToRunTimeSelectionTable(refinementController, uniformIntervals, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementControllers::uniformIntervals::uniformIntervals
(
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    refinementController(typeName, sources, dict, mesh),
    
    coeffs_(refinementDict_.optionalSubDict(typeName + "Coeffs")),
    intervals_(coeffs_.lookup<int>("intervals")),
    boundingBox_
    (
        coeffs_.lookupOrDefault<boundBox>
        (
            "boundingBox",
            boundBox(point::max, point::min)
        )
    ),
    intervalTime_(0.0),
    updateTime_(0.0)
{
    // set end time for AMR update to minimum of solution time and max beam time
    scalar endTime_ = 0.0;

    forAll(sources_, i)
    {
        endTime_ = max(sources_[i].beam().endTime(), endTime_);
    }

    intervalTime_ = min(endTime_, mesh.time().endTime().value()) / intervals_;

    Info<< "AMR time interval: " << intervalTime_ << " s" << endl;

    // TODO: this is not being used and the moment but should allow a user specified bounding box
    Info<< "refinement bounding box: " << boundingBox_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::refinementControllers::uniformIntervals::update()
{
    //- Update refinement field at update time
    if 
    (
        (mesh_.time().value() >= updateTime_)
     && refinementController::update()
    )
    {
        // update last refinement index
        lastRefinementIndex_ =
            (mesh_.time().timeIndex() == 0) ? 1 : mesh_.time().timeIndex();

        // update next refinement update time
        updateTime_ = mesh_.time().value() + intervalTime_;

        //- Set initial refinement field in base class
        //refinementController::setRefinementField();

        forAll(mesh_.cells(), celli)
        {
            refinementField_[celli] = 0.0;
        }
        refinementField_.correctBoundaryConditions();
        
        //- Calculate the bounding box for each cell
        List<treeBoundBox> cellBbs(mesh_.nCells());
        
        const pointField& points = mesh_.points();

        const vector extend = 1e-10 * vector::one;

        forAll(mesh_.cells(), celli)
        {
            treeBoundBox cellBb(point::max, point::min);

            const labelList& vertices = mesh_.cellPoints()[celli];

            forAll(vertices, j)
            {
                cellBb.min()
                    = min(cellBb.min(), points[vertices[j]] - extend);
                cellBb.max()
                    = max(cellBb.max(), points[vertices[j]] + extend);
            }

            cellBbs[celli] = cellBb;
        }

        //- Update refinement marker field:
        //-   1. March along beam paths until next update time.
        //-   2. Mark overlap cells for refinement.
        forAll(sources_, i)
        {
            const movingBeam& beam_ = sources_[i].beam();
            
            scalar time_ = mesh_.time().value();

            vector offset_ = 1.5*sources_[i].dimensions();

            while ((updateTime_ - time_) > small)
            {
                vector position_ = beam_.position(time_);

                treeBoundBox beamBb(position_ - offset_, position_ + offset_);
                
                forAll(mesh_.cells(), celli)
                {
                    if (refinementField_[celli] > 0.0)
                    {
                        // Do nothing, cell already marked for refiment
                    }
                    else if (cellBbs[celli].overlaps(beamBb))
                    {
                        refinementField_[celli] = 1.0;
                    }
                }
                refinementField_.correctBoundaryConditions();

                // stop scan path based refinement after path ends
                if ((beam_.endTime() - time_) < small)
                {
                    break;
                }

                // Calculate time step required to resolve beam motion on mesh
                label index_ = beam_.findIndex(time_);
                segment path_ = beam_.getSegment(index_);
                scalar timeToNextPath_ = path_.time() - time_;

                // if the path end time is directly hit, step to next
                while (timeToNextPath_ < small)
                {
                    index_ = index_ + 1;
                    path_ = beam_.getSegment(index_);
                    timeToNextPath_ = path_.time() - time_;
                }

                scalar dt_ = min(timeToNextPath_, max(0, updateTime_ - time_));

                if (path_.mode() == 0)
                {
                    const scalar scanTime_ =
                        sources_[i].D2sigma() / path_.parameter();

                    dt_ = min(timeToNextPath_, scanTime_);
                }

                time_ += dt_;
            }
        }
        
        return true;
    }
    else if (mesh_.time().timeIndex() < lastRefinementIndex_ + nLevels_)
    {
        //- Perform iterative refinements
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::refinementControllers::uniformIntervals::read()
{
    if (refinementController::read())
    {
        refinementDict_ = optionalSubDict(type() + "Coeffs");

        //- Mandatory entries
        refinementDict_.lookup("intervals") >> intervals_;
        refinementDict_.lookup("boundingBox") >> boundingBox_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
