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
#include "fvc.H"

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
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
    }

    return io;
}


void Foam::refinementModel::calculateCellAABBs
(
    List<treeBoundBox>& cellAABBs
) const
{
    cellAABBs.setSize(mesh_.nCells());

    const pointField& points = mesh_.points();

    const vector extend = 1e-10*vector::one;

    forAll(mesh_.cells(), celli)
    {
        treeBoundBox cellAABB(point::max, point::min);

        const labelList& vertices = mesh_.cellPoints()[celli];

        forAll(vertices, vertexi)
        {
            cellAABB.min() =
                min(cellAABB.min(), points[vertices[vertexi]] - extend);

            cellAABB.max() =
                max(cellAABB.max(), points[vertices[vertexi]] + extend);
        }

        cellAABBs[celli] = cellAABB;
    }
}


void Foam::refinementModel::scanPathFrame
(
    const point& p0,
    const point& p1,
    vector& e0,
    vector& e1,
    vector& e2
) const
{
    const vector d = p1 - p0;

    const scalar Lxy = sqrt(sqr(d.x()) + sqr(d.y()));

    if (Lxy < small)
    {
        e0 = vector(1, 0, 0);
        e1 = vector(0, 1, 0);
    }
    else
    {
        e0 = vector(-d.y()/Lxy, d.x()/Lxy, 0);

        e1 = vector( d.x()/Lxy, d.y()/Lxy, 0);
    }

    e2 = vector(0, 0, 1);
}


bool Foam::refinementModel::cellAABBOverlapsOBB
(
    const treeBoundBox& cellAABB,
    const point& centre,
    const vector& e0,
    const vector& e1,
    const vector& e2,
    const vector& L
) const
{
    const point cellCentre = 0.5*(cellAABB.min() + cellAABB.max());

    const vector cellHalfLength = 0.5*(cellAABB.max() - cellAABB.min());

    const vector d = cellCentre - centre;

    const scalar r0 =
        cellHalfLength.x()*mag(e0.x())
      + cellHalfLength.y()*mag(e0.y())
      + cellHalfLength.z()*mag(e0.z());

    if (mag(d & e0) > L.x() + r0)
    {
        return false;
    }

    const scalar r1 =
        cellHalfLength.x()*mag(e1.x())
      + cellHalfLength.y()*mag(e1.y())
      + cellHalfLength.z()*mag(e1.z());

    if (mag(d & e1) > L.y() + r1)
    {
        return false;
    }

    const scalar r2 =
        cellHalfLength.x()*mag(e2.x())
      + cellHalfLength.y()*mag(e2.y())
      + cellHalfLength.z()*mag(e2.z());

    if (mag(d & e2) > L.z() + r2)
    {
        return false;
    }

    return true;
}


Foam::scalar Foam::refinementModel::markScanPathInterval
(
    const scalar intervalStartTime,
    const scalar intervalEndTime,
    const List<treeBoundBox>& cellAABBs,
    const bool commit
)
{
    if ((intervalEndTime - intervalStartTime) <= small)
    {
        return Zero;
    }

    List<label> locallyMarked(mesh_.nCells(), 0);

    scalar addedScanPathRefinementVolume = Zero;

    forAll(sources_, sourcei)
    {
        const movingBeam& beam = sources_[sourcei].beam();

        const DynamicList<pathVector>& pathVectors = beam.path();

        if (!pathVectors.size())
        {
            continue;
        }

        const scalar projectionEndTime =
            min(intervalEndTime, beam.endTime());

        if ((projectionEndTime - intervalStartTime) <= small)
        {
            continue;
        }

        forAll(pathVectors, pathVectori)
        {
            const pathVector& pv = pathVectors[pathVectori];

            if (pv.power() <= small)
            {
                continue;
            }

            if ((pv.endTime() - intervalStartTime) <= small)
            {
                continue;
            }

            if ((projectionEndTime - pv.startTime()) <= small)
            {
                break;
            }

            const scalar t0 = max(intervalStartTime, pv.startTime());

            const scalar t1 = min(projectionEndTime, pv.endTime());

            if ((t1 - t0) <= small)
            {
                continue;
            }

            const point p0 = pv.position(t0);

            const point p1 = pv.position(t1);

            vector e0(vector::zero);
            vector e1(vector::zero);
            vector e2(vector::zero);

            scanPathFrame(p0, p1, e0, e1, e2);

            const vector d = p1 - p0;

            //- XY scan length, consistent with e1 being defined in the scan
            //  plane. 3D scan vectors are not currently supported.
            const scalar scanLength =
                sqrt(sqr(d.x()) + sqr(d.y()));

            const point centre = 0.5*(p0 + p1);

            //- OBB half-lengths.
            //  buffer_.x(): half-span along e0, width/lateral direction
            //  buffer_.y(): endpoint padding along e1, scan direction
            //  buffer_.z(): half-span along e2, depth/vertical direction
            const vector L
            (
                buffer_.x(),
                0.5*scanLength + buffer_.y(),
                buffer_.z()
            );

            const vector pathAABBHalfLength
            (
                mag(e0.x())*L.x()
              + mag(e1.x())*L.y()
              + mag(e2.x())*L.z(),

                mag(e0.y())*L.x()
              + mag(e1.y())*L.y()
              + mag(e2.y())*L.z(),

                mag(e0.z())*L.x()
              + mag(e1.z())*L.y()
              + mag(e2.z())*L.z()
            );

            treeBoundBox pathAABB
            (
                centre - pathAABBHalfLength,
                centre + pathAABBHalfLength
            );

            forAll(mesh_.cells(), celli)
            {
                if (refinementField_[celli] > 0 || locallyMarked[celli])
                {
                    continue;
                }

                if
                (
                    cellAABBs[celli].overlaps(pathAABB)
                 && cellAABBOverlapsOBB
                    (
                        cellAABBs[celli],
                        centre,
                        e0,
                        e1,
                        e2,
                        L
                    )
                )
                {
                    locallyMarked[celli] = 1;

                    addedScanPathRefinementVolume += mesh_.V()[celli];

                    if (commit)
                    {
                        refinementField_[celli] = 1;
                    }
                }
            }
            
            refinementField_.correctBoundaryConditions();
        }
    }

    reduce(addedScanPathRefinementVolume, sumOp<scalar>());

    return addedScanPathRefinementVolume;
}


Foam::scalar Foam::refinementModel::nextPoweredPathEventTime
(
    const scalar time
) const
{
    scalar nextEventTime = endTime_;

    forAll(sources_, sourcei)
    {
        const movingBeam& beam = sources_[sourcei].beam();

        const DynamicList<pathVector>& pathVectors = beam.path();

        forAll(pathVectors, pathVectori)
        {
            const pathVector& pv = pathVectors[pathVectori];

            if (pv.power() <= small)
            {
                continue;
            }

            if ((pv.startTime() - time) > small)
            {
                nextEventTime = min(nextEventTime, pv.startTime());
                break;
            }

            if ((pv.endTime() - time) > small)
            {
                nextEventTime = min(nextEventTime, pv.endTime());
                break;
            }
        }
    }

    return nextEventTime;
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
    refinementDict_(heatSourceDict_.optionalSubDict("refinementModel")),
    refine_
    (
        refinementDict_.lookupOrDefault<bool>
        (
            "refine",
            false
        )
    ),
    nLevels_
    (
        refinementDict_.lookupOrDefault<label>
        (
            "nLevels",
            0
        )
    ),
    refinementTemperature_
    (
        refinementDict_.lookupOrDefault<scalar>
        (
            "refinementTemperature",
            GREAT
        ) 
    ),
    buffer_
    (   
        refinementDict_.lookupOrDefault<vector>
        (
            "buffer",
            vector::zero
        )
    ),
    endTime_(0.0),
    refinementField_
    (
        IOobject
        (
            "refinementField",
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    )
{
    forAll(sources_, sourcei)
    {
        endTime_ = max(sources_[sourcei].beam().endTime(), endTime_);
    }

    endTime_ = min(endTime_, mesh.time().endTime().value());

    Info << "refinementModel: Performing refinement until "
         << endTime_ << " s of simulation time." << endl;

    Info << "refinementModel: Refinement buffer " << buffer_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::refinementModel::refineUsingTemperature()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    refinementField_ =
        pos0(T - dimensionedScalar(dimTemperature, refinementTemperature_));

    refinementField_.correctBoundaryConditions();
}


void Foam::refinementModel::refineUsingTime
(
    const Foam::scalar& refinementTime
)
{
    List<treeBoundBox> cellAABBs;

    calculateCellAABBs(cellAABBs);

    const scalar projectionStartTime = mesh_.time().value();

    const scalar projectionEndTime = min(endTime_, refinementTime);

    if ((projectionEndTime - projectionStartTime) > small)
    {
        markScanPathInterval
        (
            projectionStartTime,
            projectionEndTime,
            cellAABBs,
            true
        );
    }

    refinementField_.correctBoundaryConditions();
}


Foam::dimensionedScalar Foam::refinementModel::refineUsingVolume
(
    const Foam::dimensionedScalar& refinementVolume,
    const Foam::scalar&
)
{
    const scalar projectionStartTime = mesh_.time().value();

    const label maximumVolumeSearchIterations =
        refinementDict_.lookupOrDefault<label>("volumeSearchMaxIter", 10);

    const scalar volumeSearchTimeTolerance =
        refinementDict_.lookupOrDefault<scalar>
        (
            "volumeSearchTimeTolerance",
            1e-5
        );

    const scalar minRefinementVolumeFactor =
        refinementDict_.lookupOrDefault<scalar>
        (
            "minRefinementVolumeFactor",
            1.0
        );

    //- Characteristic symmetric spot volume:
    //  (2*buffer.x) * (2*buffer.y) * (2*buffer.z)
    const scalar characteristicRefinementVolume =
        8.0*buffer_.x()*buffer_.y()*buffer_.z();

    const scalar minimumScanPathRefinementVolume =
        minRefinementVolumeFactor*characteristicRefinementVolume;

    List<treeBoundBox> cellAABBs;

    calculateCellAABBs(cellAABBs);

    const scalar existingRefinementVolume =
        fvc::domainIntegrate(refinementField_).value();

    scalar refinementEndTime = projectionStartTime;

    scalar committedScanPathRefinementVolume = Zero;

    const scalar targetGlobalRefinementVolume = refinementVolume.value();

    if
    (
        existingRefinementVolume >= targetGlobalRefinementVolume
     && minimumScanPathRefinementVolume <= small
    )
    {
        Info << "Existing global refinement volume already satisfies target."
             << endl;

        Info << "Minimum scan-path refinement volume: "
             << minimumScanPathRefinementVolume << endl;

        Info << "Actual scan-path refinement volume: "
             << committedScanPathRefinementVolume << endl;

        Info << "Actual global refinement volume: "
             << existingRefinementVolume << endl;

        return dimensionedScalar(dimTime, refinementEndTime);
    }

    while ((endTime_ - refinementEndTime) > small)
    {
        scalar nextProjectionEndTime =
            min(nextPoweredPathEventTime(refinementEndTime), endTime_);

        if ((nextProjectionEndTime - refinementEndTime) <= small)
        {
            break;
        }

        const scalar trialScanPathRefinementVolume =
            markScanPathInterval
            (
                refinementEndTime,
                nextProjectionEndTime,
                cellAABBs,
                false
            );

        const scalar trialTotalScanPathRefinementVolume =
            committedScanPathRefinementVolume
          + trialScanPathRefinementVolume;

        const scalar trialGlobalRefinementVolume =
            existingRefinementVolume
          + trialTotalScanPathRefinementVolume;

        const bool targetReached =
        (
            trialGlobalRefinementVolume >= targetGlobalRefinementVolume
         && trialTotalScanPathRefinementVolume
            >= minimumScanPathRefinementVolume
        );

        if (targetReached)
        {
            scalar lowerRefinementTime = refinementEndTime;

            scalar upperRefinementTime = nextProjectionEndTime;

            for
            (
                label iteration = 0;
                iteration < maximumVolumeSearchIterations;
                ++iteration
            )
            {
                if
                (
                    (upperRefinementTime - lowerRefinementTime)
                 <= volumeSearchTimeTolerance
                )
                {
                    break;
                }

                const scalar midpointRefinementTime =
                    0.5*(lowerRefinementTime + upperRefinementTime);

                const scalar midpointScanPathRefinementVolume =
                    markScanPathInterval
                    (
                        refinementEndTime,
                        midpointRefinementTime,
                        cellAABBs,
                        false
                    );

                const scalar midpointTotalScanPathRefinementVolume =
                    committedScanPathRefinementVolume
                  + midpointScanPathRefinementVolume;

                const scalar midpointGlobalRefinementVolume =
                    existingRefinementVolume
                  + midpointTotalScanPathRefinementVolume;

                if
                (
                    midpointGlobalRefinementVolume
                        >= targetGlobalRefinementVolume
                 && midpointTotalScanPathRefinementVolume
                        >= minimumScanPathRefinementVolume
                )
                {
                    upperRefinementTime = midpointRefinementTime;
                }
                else
                {
                    lowerRefinementTime = midpointRefinementTime;
                }
            }

            committedScanPathRefinementVolume +=
                markScanPathInterval
                (
                    refinementEndTime,
                    upperRefinementTime,
                    cellAABBs,
                    true
                );

            refinementEndTime = upperRefinementTime;

            break;
        }

        committedScanPathRefinementVolume +=
            markScanPathInterval
            (
                refinementEndTime,
                nextProjectionEndTime,
                cellAABBs,
                true
            );

        refinementEndTime = nextProjectionEndTime;
    }

    refinementField_.correctBoundaryConditions();

    Info << "Minimum scan-path refinement volume: "
         << minimumScanPathRefinementVolume << endl;

    Info << "Actual scan-path refinement volume: "
         << committedScanPathRefinementVolume << endl;

    Info << "Actual global refinement volume: "
         << existingRefinementVolume + committedScanPathRefinementVolume
         << endl;

    Info << "Refinement volume search ended at time: "
         << refinementEndTime << endl;

    return dimensionedScalar(dimTime, refinementEndTime);
}


bool Foam::refinementModel::read()
{
    return regIOobject::read();
}

// ************************************************************************* //
