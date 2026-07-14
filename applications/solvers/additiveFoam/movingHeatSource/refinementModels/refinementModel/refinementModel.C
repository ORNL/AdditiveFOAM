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
        mesh,
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

Foam::label Foam::refinementModel::readMaxRefinementLevel() const
{
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    return dynamicMeshDict.subDict("topoChanger").lookup<label>
    (
        "maxRefinement"
    );
}


void Foam::refinementModel::readSourceBuffers()
{
    sourceBuffers_.setSize(sources_.size());

    const dictionary& buffersDict = refinementDict_.subDict("buffers");

    forAll(sources_, sourcei)
    {
        buffersDict.lookup(sources_[sourcei].sourceName())
            >> sourceBuffers_[sourcei];
    }
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
        //- Scan direction
        e0 = vector(d.x()/Lxy, d.y()/Lxy, 0);

        //- Transverse direction
        e1 = vector(-d.y()/Lxy, d.x()/Lxy, 0);
    }

    //- Depth/vertical direction
    e2 = vector(0, 0, 1);
}


bool Foam::refinementModel::cellOverlapsOBB
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


Foam::scalar Foam::refinementModel::markScanPath
(
    const scalar startTime,
    const scalar endTime,
    const List<treeBoundBox>& cellAABBs,
    const bool commit
)
{
    if ((endTime - startTime) <= small)
    {
        return Zero;
    }

    List<label> locallyMarked(mesh_.nCells(), 0);

    scalar addedScanPathRefineVolume = Zero;

    forAll(sources_, sourcei)
    {
        const movingBeam& beam = sources_[sourcei].beam();

        const DynamicList<pathVector>& pathVectors = beam.path();

        if (!pathVectors.size())
        {
            continue;
        }

        const scalar projectedEndTime = min(endTime, beam.endTime());

        if ((projectedEndTime - startTime) <= small)
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

            if ((pv.endTime() - startTime) <= small)
            {
                continue;
            }

            if ((projectedEndTime - pv.startTime()) <= small)
            {
                break;
            }

            const scalar t0 = max(startTime, pv.startTime());

            const scalar t1 = min(projectedEndTime, pv.endTime());

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

            const scalar scanLength = sqrt(sqr(d.x()) + sqr(d.y()));

            const point centre = 0.5*(p0 + p1);

            const vector L
            (
                0.5*scanLength + sourceBuffers_[sourcei].x(),
                sourceBuffers_[sourcei].y(),
                sourceBuffers_[sourcei].z()
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
                 && cellOverlapsOBB
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

                    addedScanPathRefineVolume += mesh_.V()[celli];

                    if (commit)
                    {
                        refinementField_[celli] = 1;
                    }
                }
            }
        }
    }

    reduce(addedScanPathRefineVolume, sumOp<scalar>());

    return addedScanPathRefineVolume;
}


Foam::scalar Foam::refinementModel::nextPoweredPathEventTime
(
    const scalar time
) const
{
    scalar nextEventTime = scanEndTime_;

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


Foam::scalar Foam::refinementModel::sourceBufferVolume() const
{
    scalar bufferVolume = Zero;

    forAll(sourceBuffers_, sourcei)
    {
        bufferVolume += 8.0*cmptProduct(sourceBuffers_[sourcei]);
    }

    return bufferVolume;
}


Foam::scalar Foam::refinementModel::minimumScanPathRefineVolume
(
    const scalar nBufferVolumes
) const
{
    return nBufferVolumes*sourceBufferVolume();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementModel::refinementModel
(
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
    maxRefinementLevel_(0),
    refinementTemperature_(GREAT),
    sourceBuffers_(sources.size(), vector::zero),
    scanEndTime_(0.0),
    refinementField_
    (
        IOobject
        (
            "refinementField",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    )
{}


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
    refinementDict_(heatSourceDict_.subDict("refinementModel")),
    maxRefinementLevel_(readMaxRefinementLevel()),
    refinementTemperature_
    (
        refinementDict_.lookupOrDefault<scalar>("refinementTemperature", GREAT)
    ),
    sourceBuffers_(sources.size(), vector::zero),
    scanEndTime_(0.0),
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
        scanEndTime_ = max(sources_[sourcei].beam().endTime(), scanEndTime_);
    }

    scanEndTime_ = min(scanEndTime_, mesh.time().endTime().value());

    readSourceBuffers();

    Info << "refinementModel: Performing refinement until "
         << scanEndTime_ << " s of simulation time." << endl;

    forAll(sources_, sourcei)
    {
        Info << "refinementModel: Refinement buffer for "
             << sources_[sourcei].sourceName() << " "
             << sourceBuffers_[sourcei] << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::refinementModel::markTemperature()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    refinementField_ =
        pos0(T - dimensionedScalar(dimTemperature, refinementTemperature_));

    refinementField_.correctBoundaryConditions();
}


void Foam::refinementModel::markScanPathTime
(
    const Foam::scalar& refineTime
)
{
    List<treeBoundBox> cellAABBs;

    calculateCellAABBs(cellAABBs);

    const scalar projectionStartTime = mesh_.time().value();

    const scalar projectedEndTime = min(scanEndTime_, refineTime);

    if ((projectedEndTime - projectionStartTime) > small)
    {
        markScanPath
        (
            projectionStartTime,
            projectedEndTime,
            cellAABBs,
            true
        );
    }

    refinementField_.correctBoundaryConditions();
}


Foam::dimensionedScalar Foam::refinementModel::markScanPathVolume
(
    const Foam::dimensionedScalar& targetRefineVolume,
    const scalar minScanPathRefineVolume,
    const label maxSearchIter,
    const scalar timeTolerance
)
{
    const scalar projectionStartTime = mesh_.time().value();

    List<treeBoundBox> cellAABBs;

    calculateCellAABBs(cellAABBs);

    const scalar existingRefineVolume =
        fvc::domainIntegrate(refinementField_).value();

    scalar refineEndTime = projectionStartTime;

    scalar committedScanPathRefineVolume = Zero;

    const scalar targetGlobalRefineVolume = targetRefineVolume.value();

    const bool globalTargetReached =
    (
        (existingRefineVolume >= targetGlobalRefineVolume)
     && (minScanPathRefineVolume <= small)
    );

    if (globalTargetReached)
    {
        Info << "Existing global refine volume already satisfies target."
             << endl;

        Info << "Min scan-path refine volume: "
             << minScanPathRefineVolume << endl;

        Info << "Actual scan-path refine volume: "
             << committedScanPathRefineVolume << endl;

        Info << "Actual global refine volume: "
             << existingRefineVolume << endl;

        return dimensionedScalar(dimTime, refineEndTime);
    }

    while ((scanEndTime_ - refineEndTime) > small)
    {
        scalar nextProjectedEndTime =
            min(nextPoweredPathEventTime(refineEndTime), scanEndTime_);

        if ((nextProjectedEndTime - refineEndTime) <= small)
        {
            break;
        }

        const scalar trialScanPathRefineVolume =
            markScanPath
            (
                refineEndTime,
                nextProjectedEndTime,
                cellAABBs,
                false
            );

        const scalar trialTotalScanPathRefineVolume =
            committedScanPathRefineVolume
          + trialScanPathRefineVolume;

        const scalar trialGlobalRefineVolume =
            existingRefineVolume
          + trialTotalScanPathRefineVolume;

        const bool targetReached =
        (
            (trialGlobalRefineVolume >= targetGlobalRefineVolume)
         && (trialTotalScanPathRefineVolume >= minScanPathRefineVolume)
        );

        if (targetReached)
        {
            scalar lowerRefineTime = refineEndTime;

            scalar upperRefineTime = nextProjectedEndTime;

            for
            (
                label iteration = 0;
                iteration < maxSearchIter;
                ++iteration
            )
            {
                if
                (
                    (upperRefineTime - lowerRefineTime)
                 <= timeTolerance
                )
                {
                    break;
                }

                const scalar midpointRefineTime =
                    0.5*(lowerRefineTime + upperRefineTime);

                const scalar midpointScanPathRefineVolume =
                    markScanPath
                    (
                        refineEndTime,
                        midpointRefineTime,
                        cellAABBs,
                        false
                    );

                const scalar midpointTotalScanPathRefineVolume =
                    committedScanPathRefineVolume
                  + midpointScanPathRefineVolume;

                const scalar midpointGlobalRefineVolume =
                    existingRefineVolume
                  + midpointTotalScanPathRefineVolume;

                const bool midpointTargetReached =
                (
                    (
                        midpointGlobalRefineVolume
                     >= targetGlobalRefineVolume
                    )
                 && (
                        midpointTotalScanPathRefineVolume
                     >= minScanPathRefineVolume
                    )
                );

                Info<< "refinementModel: volume search iteration "
                    << iteration << nl
                    << "    lowerRefineTime: " << lowerRefineTime << nl
                    << "    midpointRefineTime: " << midpointRefineTime << nl
                    << "    upperRefineTime: " << upperRefineTime << nl
                    << "    midpoint scan-path refine volume: "
                    << midpointTotalScanPathRefineVolume << nl
                    << "    midpoint global refine volume: "
                    << midpointGlobalRefineVolume << nl
                    << "    target reached: "
                    << midpointTargetReached << endl;

                if (midpointTargetReached)
                {
                    upperRefineTime = midpointRefineTime;
                }
                else
                {
                    lowerRefineTime = midpointRefineTime;
                }
            }

            committedScanPathRefineVolume +=
                markScanPath
                (
                    refineEndTime,
                    upperRefineTime,
                    cellAABBs,
                    true
                );

            refineEndTime = upperRefineTime;

            break;
        }

        committedScanPathRefineVolume +=
            markScanPath
            (
                refineEndTime,
                nextProjectedEndTime,
                cellAABBs,
                true
            );

        refineEndTime = nextProjectedEndTime;
    }

    refinementField_.correctBoundaryConditions();

    Info << "Min scan-path refine volume: "
         << minScanPathRefineVolume << endl;

    Info << "Actual scan-path refine volume: "
         << committedScanPathRefineVolume << endl;

    Info << "Actual global refine volume: "
         << existingRefineVolume + committedScanPathRefineVolume << endl;

    Info << "Refine volume search ended at time: "
         << refineEndTime << endl;

    return dimensionedScalar(dimTime, refineEndTime);
}


bool Foam::refinementModel::read()
{
    if (regIOobject::read())
    {
        refinementDict_ = optionalSubDict("refinementModel");

        maxRefinementLevel_ = readMaxRefinementLevel();

        refinementTemperature_ =
            refinementDict_.lookupOrDefault<scalar>
            (
                "refinementTemperature",
                GREAT
            );

        readSourceBuffers();

        return true;
    }

    return false;
}

// ************************************************************************* //
