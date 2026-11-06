/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ExaCA.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "pointMVCWeight.H"
#include "Pstream.H"

#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "thermoPath.H"

#include <cmath>

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ExaCA, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ExaCA,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ExaCA::ExaCA
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    T_(mesh_.lookupObject<volScalarField>("T")),
    meshChanged_(true),
    nPoints_(labelVector(vector::one)),
    box_(point::max, point::min),
    isoValue_(Zero),
    dx_(Zero)
{
    read(dict);

    updateMeshState();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ExaCA::~ExaCA()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::functionObjects::ExaCA::resetMeshState()
{
    cellMapIds.clear();

    overlapCells.clear();

    meshChanged_ = true;
}


void Foam::functionObjects::ExaCA::updateMeshState()
{
    cellMapIds.clear();

    overlapCells.clear();

    setOverlapCells();

    meshChanged_ = false;
}


Foam::label Foam::functionObjects::ExaCA::pointId
(
    const label i,
    const label j,
    const label k
) const
{
    return i + nPoints_.x()*(j + nPoints_.y()*k);
}


Foam::label Foam::functionObjects::ExaCA::getCellMap(const label celli)
{
    if (cellMapIds.found(celli))
    {
        return cellMapIds[celli];
    }

    const label mapId = cellMaps.size();

    cellMapIds.set(celli, mapId);

    cellMaps.append(makeCellMap(celli));

    return mapId;
}


Foam::functionObjects::ExaCA::cellMap
Foam::functionObjects::ExaCA::makeCellMap(const label celli) const
{
    cellMap cm;

    const pointField& points = mesh_.points();

    const labelList& vertices = mesh_.cellPoints()[celli];

    boundBox cellBb(point::max, point::min);

    const vector extend = 1e-10*vector::one;

    forAll(vertices, vertexi)
    {
        const point& p = points[vertices[vertexi]];

        cellBb.min() = min(cellBb.min(), p - extend);
        cellBb.max() = max(cellBb.max(), p + extend);
    }

    if (!cellBb.overlaps(box_))
    {
        return cm;
    }

    const point& boxMax = box_.max();

    const label i0 =
        max
        (
            label(0),
            label(std::ceil((boxMax.x() - cellBb.max().x())/dx_))
        );

    const label i1 =
        min
        (
            nPoints_.x() - 1,
            label(std::floor((boxMax.x() - cellBb.min().x())/dx_))
        );

    const label j0 =
        max
        (
            label(0),
            label(std::ceil((boxMax.y() - cellBb.max().y())/dx_))
        );

    const label j1 =
        min
        (
            nPoints_.y() - 1,
            label(std::floor((boxMax.y() - cellBb.min().y())/dx_))
        );

    const label k0 =
        max
        (
            label(0),
            label(std::ceil((boxMax.z() - cellBb.max().z())/dx_))
        );

    const label k1 =
        min
        (
            nPoints_.z() - 1,
            label(std::floor((boxMax.z() - cellBb.min().z())/dx_))
        );

    if (i1 < i0 || j1 < j0 || k1 < k0)
    {
        return cm;
    }

    const cell& cFaces = mesh_.cells()[celli];

    const label nPlanes = cFaces.size();

    scalarField nx(nPlanes);
    scalarField ny(nPlanes);
    scalarField nz(nPlanes);
    scalarField d(nPlanes);

    const vectorField& faceAreas = mesh_.faceAreas();
    const pointField& faceCentres = mesh_.faceCentres();
    const labelUList& faceOwner = mesh_.faceOwner();

    forAll(cFaces, cFacei)
    {
        const label facei = cFaces[cFacei];

        const vector& Sf = faceAreas[facei];
        const point& Cf = faceCentres[facei];

        if (faceOwner[facei] == celli)
        {
            nx[cFacei] = Sf.x();
            ny[cFacei] = Sf.y();
            nz[cFacei] = Sf.z();
        }
        else
        {
            nx[cFacei] = -Sf.x();
            ny[cFacei] = -Sf.y();
            nz[cFacei] = -Sf.z();
        }

        d[cFacei] =
            nx[cFacei]*Cf.x()
          + ny[cFacei]*Cf.y()
          + nz[cFacei]*Cf.z();
    }

    const scalar xMax = boxMax.x();
    const scalar yMax = boxMax.y();
    const scalar zMax = boxMax.z();

    const scalar eps = 1e-10;

    for (label k = k0; k <= k1; ++k)
    {
        const scalar z = zMax - scalar(k)*dx_;
        const scalar sz = z - eps;

        for (label j = j0; j <= j1; ++j)
        {
            const scalar y = yMax - scalar(j)*dx_;
            const scalar sy = y - eps;

            for (label i = i0; i <= i1; ++i)
            {
                const scalar x = xMax - scalar(i)*dx_;
                const scalar sx = x - eps;

                bool pointInCell = true;

                for (label planei = 0; planei < nPlanes; ++planei)
                {
                    if
                    (
                        nx[planei]*sx
                      + ny[planei]*sy
                      + nz[planei]*sz
                      - d[planei]
                      > 0
                    )
                    {
                        pointInCell = false;
                        break;
                    }
                }

                if (!pointInCell)
                {
                    continue;
                }

                const point pt(x, y, z);

                const point spt(sx, sy, sz);

                pointMVCWeight cpw(mesh_, spt, celli);

                cm.pointIds.append(pointId(i, j, k));

                cm.positions.append(pt);

                cm.weights.append(cpw.weights());
            }
        }
    }

    cm.pointIds.shrink();

    cm.positions.shrink();

    cm.weights.shrink();

    return cm;
}


void Foam::functionObjects::ExaCA::appendInterval
(
    const label celli,
    const pointScalarField& Tp0,
    const pointScalarField& Tp,
    const scalar t0,
    const scalar t
)
{
    const label mapId = getCellMap(celli);

    if (!cellMaps[mapId].pointIds.size())
    {
        return;
    }

    const labelList& vertices = mesh_.cellPoints()[celli];

    eventInterval event;

    event.mapId = mapId;

    event.t0 = t0;

    event.t1 = t;

    event.psi0.setSize(vertices.size());

    event.psi.setSize(vertices.size());

    forAll(vertices, pI)
    {
        event.psi0[pI] = Tp0[vertices[pI]];

        event.psi[pI] = Tp[vertices[pI]];
    }

    intervals.append(event);
}


void Foam::functionObjects::ExaCA::writeData() const
{
    const fileName exacaPath
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()/"ExaCA"
    );

    mkDir(exacaPath);

    OFstream os
    (
        exacaPath + "/" + "data_" + Foam::name(Pstream::myProcNo()) + ".csv"
    );

    os << "x,y,z,tm,ts,cr" << endl;

    forAll(data, i)
    {
        const List<scalar>& row = data[i];

        os  << row[0] << ","
            << row[1] << ","
            << row[2] << ","
            << row[3] << ","
            << row[4] << ","
            << row[5] << "\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ExaCA::read(const dictionary& dict)
{
    box_ = dict.lookup("box");

    if (!dict.readIfPresent("isoValue", isoValue_))
    {
        isoValue_ = thermoPath(mesh_).liquidus();
    }

    dx_ = dict.lookup<scalar>("dx");

    nPoints_ = labelVector(vector::one + box_.span()/dx_);

    return true;
}


void Foam::functionObjects::ExaCA::setOverlapCells()
{
    overlapCells.clear();

    const pointField& points = mesh_.points();

    boundBox procBb(points);

    const vector extend = 1e-10*vector::one;

    if (procBb.overlaps(box_))
    {
        forAll(mesh_.cells(), celli)
        {
            boundBox cellBb(point::max, point::min);

            const labelList& vertices = mesh_.cellPoints()[celli];

            forAll(vertices, i)
            {
                cellBb.min() = min(cellBb.min(), points[vertices[i]] - extend);
                cellBb.max() = max(cellBb.max(), points[vertices[i]] + extend);
            }

            if (cellBb.overlaps(box_))
            {
                overlapCells.append(celli);
            }
        }
    }

    overlapCells.shrink();
}


Foam::wordList Foam::functionObjects::ExaCA::fields() const
{
    return wordList(1, T_.name());
}


bool Foam::functionObjects::ExaCA::execute()
{
    if (meshChanged_)
    {
        updateMeshState();
    }

    const volPointInterpolation& vpi = volPointInterpolation::New(mesh_);

    const pointScalarField Tp0
    (
        IOobject
        (
            "Tp0_",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vpi.interpolate(T_.oldTime())
    );

    const pointScalarField Tp
    (
        IOobject
        (
            "Tp_",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vpi.interpolate(T_)
    );

    const scalar t = mesh_.time().value();

    const scalar t0 = t - mesh_.time().deltaTValue();

    forAll(overlapCells, i)
    {
        const label celli = overlapCells[i];

        const labelList& vertices = mesh_.cellPoints()[celli];

        label c0 = 0;

        label c1 = 0;

        forAll(vertices, pI)
        {
            if (Tp0[vertices[pI]] >= isoValue_)
            {
                c0++;
            }

            if (Tp[vertices[pI]] >= isoValue_)
            {
                c1++;
            }
        }

        const label n = vertices.size();

        if (!(c0 % n) && !(c1 % n))
        {
            if (c0 != c1)
            {
                c0 = 0;

                c1 = 1;
            }
        }

        c0 %= n;

        c1 %= n;

        if (c0 || c1)
        {
            appendInterval(celli, Tp0, Tp, t0, t);
        }
    }

    return true;
}


void Foam::functionObjects::ExaCA::interpolate()
{
    if (intervals.size() == 0)
    {
        writeData();

        return;
    }

    meltTimes.clear();

    data.clear();

    forAll(intervals, eventi)
    {
        const eventInterval& event = intervals[eventi];

        if (event.t1 <= event.t0 + small)
        {
            continue;
        }

        const cellMap& cm = cellMaps[event.mapId];

        forAll(cm.pointIds, pointi)
        {
            scalar tp0 = Zero;

            scalar tp = Zero;

            const scalarField& w = cm.weights[pointi];

            forAll(w, j)
            {
                tp0 += w[j]*event.psi0[j];

                tp += w[j]*event.psi[j];
            }

            const label caPointId = cm.pointIds[pointi];

            if ((tp <= isoValue_) && (tp0 > isoValue_))
            {
                const point& pt = cm.positions[pointi];

                const scalar m =
                    min
                    (
                        max
                        (
                            (isoValue_ - tp0)/(tp - tp0),
                            scalar(0)
                        ),
                        scalar(1)
                    );

                scalar tm = event.t0;

                if (meltTimes.found(caPointId))
                {
                    tm = meltTimes[caPointId];
                }

                data.append
                (
                    {
                        pt[0],
                        pt[1],
                        pt[2],
                        tm,
                        event.t0 + m*(event.t1 - event.t0),
                        (tp0 - tp)/(event.t1 - event.t0)
                    }
                );
            }
            else if ((tp > isoValue_) && (tp0 <= isoValue_))
            {
                const scalar m =
                    min
                    (
                        max
                        (
                            (isoValue_ - tp0)/(tp - tp0),
                            scalar(0)
                        ),
                        scalar(1)
                    );

                meltTimes.set(caPointId, event.t0 + m*(event.t1 - event.t0));
            }
        }
    }

    data.shrink();

    writeData();
}


bool Foam::functionObjects::ExaCA::end()
{
    intervals.shrink();

    cellMaps.shrink();

    Info<< "Number of ExaCA event intervals: "
        << returnReduce(intervals.size(), sumOp<label>()) << endl;

    Info<< "Number of cached ExaCA cell maps: "
        << returnReduce(cellMaps.size(), sumOp<label>()) << endl;

    mesh_.time().cpuTimeIncrement();

    interpolate();

    Info<< "Number of solidification events: "
        << returnReduce(data.size(), sumOp<label>()) << endl;

    Info<< "Successfully interpolated and wrote ExaCA data in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>())
        << " s" << endl << endl;

    return true;
}


bool Foam::functionObjects::ExaCA::write()
{
    return true;
}


void Foam::functionObjects::ExaCA::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        resetMeshState();
    }
}


void Foam::functionObjects::ExaCA::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        resetMeshState();
    }
}


void Foam::functionObjects::ExaCA::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        resetMeshState();
    }
}


void Foam::functionObjects::ExaCA::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        FatalErrorInFunction
            << "ExaCA does not currently support dynamic mesh redistribution. "
            << "Moving meshes and topology-changing meshes on a fixed "
            << "decomposition are supported."
            << abort(FatalError);
    }
}


// ************************************************************************* //
