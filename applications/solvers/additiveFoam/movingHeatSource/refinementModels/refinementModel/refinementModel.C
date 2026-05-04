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

#include "refinementModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(refinementModel, 0);
    defineRunTimeSelectionTable(refinementModel, dictionary);
}

const Foam::word Foam::refinementModel::refinementModelDictName
(
    "refinementModelDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::refinementModel::createIOobject
(
    const dictionary& dict,
    const fvMesh& mesh
) const
{
    typeIOobject<IOdictionary> io
    (
        dict.name(),
        mesh.time().constant(),
        mesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModel::refinementModel
(
    const word& type,
    const PtrList<heatSourceModel>& sources,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(dict, mesh)),

    sources_(sources),
    mesh_(mesh),
    heatSourceDict_(dict),
    refinementDict_(heatSourceDict_.optionalSubDict("refinementControl")),
    refine_
    (
        (type != "none")
      ? refinementDict_.lookup<bool>("refine")
      : false
    ),
    nLevels_
    (
        (type != "none")
      ? refinementDict_.lookup<label>("nLevels")
      : 0
    ),
    refinementTemperature_
    (
        (type != "none")
      ? refinementDict_.lookupOrDefault<scalar>("refinementTemperature", GREAT)
      : GREAT
    ),
    buffer_
    (
        (type != "none")
        ? refinementDict_.lookupOrDefault<vector>("buffer", vector::zero)
        : vector::zero
    ),
    endTime_(0.0),
    refinementField_
    (
        IOobject
        (
            "refinementField",
            mesh_.time().name(),
            mesh_,
            refine_ ? IOobject::READ_IF_PRESENT : IOobject::NO_READ,
            refine_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    )
{
    //- Set AMR update end time to minimum of solution time and max beam time
    forAll(sources_, i)
    {
        const movingBeam& beam_ = sources_[i].beam();
        
        endTime_ = max(beam_.endTime(), endTime_);
    }

    endTime_ = min(endTime_, mesh.time().endTime().value());

    Info << "refinementModel: Performing refinement until "
         << endTime_ << " s of simulation time." << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::refinementModel::refineUsingTemperature()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    const dimensionedScalar Tr(dimTemperature, refinementTemperature_);

    refinementField_ = pos0(T - Tr);

    refinementField_.correctBoundaryConditions();
}

void Foam::refinementModel::refineUsingTime(const Foam::scalar& refineTime)
{
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
    
    //- Update the refinement marker field
    forAll(sources_, i)
    {
        const movingBeam& beam_ = sources_[i].beam();
        
        scalar time_ = mesh_.time().value();

        vector offset_ = max(buffer_, 1.5*sources_[i].dimensions());

        while ((min(beam_.endTime(), refineTime) - time_) > small)
        {
            vector position_ = beam_.position(time_);

            treeBoundBox beamBb
            (
                position_ - offset_,
                position_ + offset_
            );
            
            forAll(mesh_.cells(), celli)
            {
                if (refinementField_[celli] > 0)
                {
                    // Do nothing, cell already marked for refiment
                }
                else if (cellBbs[celli].overlaps(beamBb))
                {
                    refinementField_[celli] = 1;
                }
            }
            
            refinementField_.correctBoundaryConditions();

            //- Calculate time step required to resolve beam motion on mesh
            label index_ = beam_.findIndex(time_);
            segment path_ = beam_.getSegment(index_);
            scalar timeToNextPath_ = path_.time() - time_;

            //- If the path end time is directly hit, step to next path
            while (mag(timeToNextPath_) < small)
            {
                index_ = index_ + 1;
                path_ = beam_.getSegment(index_);
                timeToNextPath_ = path_.time() - time_;
            }

            scalar dt_ = min(timeToNextPath_, max(0, refineTime - time_));

            if (path_.mode() == 0)
            {
                const scalar scanTime_ =
                    sources_[i].D2sigma() / path_.parameter();

                dt_ = min(timeToNextPath_, scanTime_);
            }

            time_ += dt_;
        }
    }
    
    return;
}


Foam::dimensionedScalar Foam::refinementModel::refineUsingVolume
(
    const Foam::dimensionedScalar& refineVol,
    const Foam::scalar& minRefVol
)
{    
    const Foam::scalar minIntervalTime = 0.0;
    
    //- Set next refinement time to current time
    scalar refTime = mesh_.time().value();
    
    //- Find minimum refinement time
    scalar minRefTime = refTime + minIntervalTime;

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

    //- Mark current positions of all beams for refinement and find
    //  minimum time step for all beams
    scalar dt = 0.0;
    scalar refVol = 0.0;

    forAll(sources_, i)
    {
        //- Mark refinement field for beam at current time
        const movingBeam& beam = sources_[i].beam();

        vector offset = max(buffer_, 1.5 * sources_[i].dimensions());

        vector position = beam.position(refTime);

        treeBoundBox beamBb
        (
            position - offset,
            position + offset
        );

        forAll(mesh_.cells(), celli)
        {
            if (refinementField_[celli] > 0)
            {
                // Do nothing, cell already marked for refiment
            }
            else if (cellBbs[celli].overlaps(beamBb))
            {
                refinementField_[celli] = 1;
		refVol += mesh_.V()[celli];
            }
        }

        refinementField_.correctBoundaryConditions();

        //- Find minimum time step
        label index = beam.findIndex(refTime);
        segment path = beam.getSegment(index);
        scalar timeToNextPath = path.time() - refTime;

        //- If the path end time is directly hit, step to next path
        while (mag(timeToNextPath) < small)
        {
            index = index + 1;
            path = beam.getSegment(index);
            timeToNextPath = path.time() - refTime;
        }

        dt = timeToNextPath;

        if (path.mode() == 0)
        {
            const scalar scanTime =
                sources_[i].D2sigma() / path.parameter();

            dt = min(dt, scanTime);
        }
    }

    // Get volume of cells refined by current position only
    reduce(refVol, sumOp<scalar>());

    Info << "Refinement volume of current position only: " << refVol << endl;

    //- Get volume of cells refined from other criteria, i.e. melt pool,
    //  apart from those which were marked from the current beam position
    const scalar meltPoolVol = fvc::domainIntegrate(refinementField_).value() - refVol;

    Info << "Refinement volume of melt pool only: " << meltPoolVol << endl;

    //- March along scan path(s) and refine until target volume is reached
    while
    (
        ((refVol + meltPoolVol) < refineVol.value())
     || (refTime < minRefTime)
     || (refVol < minRefVol)
    )
    {
        //- Check that end time not reached
        if (refTime > endTime_)
        {
            break;
        }

        scalar refVoli = 0.0;

        forAll(sources_, i)
        {
            const movingBeam& beam = sources_[i].beam();

            vector offset = max(buffer_, 1.5 * sources_[i].dimensions());

            vector position = beam.position(refTime);

            treeBoundBox beamBb
            (
                position - offset,
                position + offset
            );

            forAll(mesh_.cells(), celli)
            {
                if (refinementField_[celli] > 0)
                {
                    // Do nothing, cell already marked for refiment
                }
                else if (cellBbs[celli].overlaps(beamBb))
                {
                    refinementField_[celli] = 1;
                    refVoli += mesh_.V()[celli];
                }
            }
        }

        reduce(refVoli, sumOp<scalar>());

        refVol += refVoli;

        refTime += dt;
    }

    refinementField_.correctBoundaryConditions();

    Info << "Actual refinement volume: " << refVol << endl;

    return dimensionedScalar(dimTime, refTime);
}


bool Foam::refinementModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
