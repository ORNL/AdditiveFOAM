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

#include "foamToExaCA/foamToExaCA.H"
#include "addToRunTimeSelectionTable.H"
#include "DynamicList.H"
#include "interpolation.H"
#include "labelVector.H"
#include "pointMVCWeight.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(foamToExaCA, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamToExaCA::foamToExaCA
(
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "foamToExaCADict",
            T.time().constant(),
            T.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),

    mesh_(T.mesh()),

    runTime_(T.time()),

    T_(T),

    vpi_(mesh_),

    Tp_
    (
        IOobject
        (
            "Tp_",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vpi_.interpolate(T_)
    ),

    execute_(false)
{
    if (this->headerOk())
    {
        execute_ = this->lookup<Switch>("execute");
    }

    initialize();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::foamToExaCA::initialize()
{
    if (!execute_)
    {
        return;
    }

    box_ = this->lookup("box");
    iso_ = this->lookup<scalar>("isotherm");
    dx_  = this->lookup<scalar>("dx");

    treeBoundBox bb(box_);

    // create a compact cell-stencil using the overlap sub-space
    const pointField& points = mesh_.points();

    treeBoundBox procBb(points);

    const vector extend = 1e-10*vector::one;

    if (procBb.overlaps(bb))
    {
        forAll(mesh_.cells(), celli)
        {
            treeBoundBox cellBb(point::max, point::min);

            const labelList& vertices = mesh_.cellPoints()[celli];

            forAll(vertices, i)
            {
                cellBb.min() = min(cellBb.min(), points[vertices[i]] - extend);
                cellBb.max() = max(cellBb.max(), points[vertices[i]] + extend);
            }

            if (cellBb.overlaps(bb))
            {
                compactCells.append(celli);
            }
        }
    }

    compactCells.shrink();
}

void Foam::foamToExaCA::update()
{
    if (!execute_)
    {
        return;
    }

    const pointScalarField Tp0_("Tp0_", Tp_);

    Tp_ = vpi_.interpolate(T_);

    // capture events at interface cells: {cell id, time, vertice temperatures}
    forAll(compactCells, i)
    {
        label celli = compactCells[i];

        const labelList& vertices = mesh_.cellPoints()[celli];

        label prev = 0;
        label curr = 0;

        forAll(vertices, pI)
        {
            if (Tp0_[vertices[pI]] >= iso_)
            {
                prev++;
            }

            if (Tp_[vertices[pI]] >= iso_)
            {
                curr++;
            }
        }

        prev = prev % vertices.size();
        curr = curr % vertices.size();

        // overshoot correction: interface jumped this cell during time step
        if (!prev && !curr)
        {
            if
            (
                (T_[celli] > iso_ && T_.oldTime()[celli] < iso_)
             || (T_[celli] < iso_ && T_.oldTime()[celli] > iso_)
            )
            {
                prev = 0;
                curr = 1;
            }
        }

        // capture the solidification events
        if (prev || curr)
        {
            List<scalar> event(vertices.size() + 2);

            event[0] = celli;
        
            // add previous event
            if (curr && !prev)
            {
                event[1] = runTime_.value() - runTime_.deltaTValue();

                forAll(vertices, pI)
                {
                    event[pI + 2] = Tp0_[vertices[pI]];
                }

                events.append(event);
            }

            // add current event
            event[1] = runTime_.value();

            forAll(vertices, pI)
            {
                event[pI + 2] = Tp_[vertices[pI]];
            }

            events.append(event);
        }
    }
}

void Foam::foamToExaCA::interpolatePoints()
{
    // dynamic list for exaca reduced data format
    DynamicList<List<scalar>> data;

    // initialize melting time for first event
    List<scalar> tm;
    tm.setSize(pointsInCell[events[0][0]].size(), events[0][1]);

    // events are order by cell id, and each cell id is ordered in time
    for (label i = 1; i < events.size(); i++)
    {
        const List<scalar> prevEvent = events[i - 1];

        const List<scalar> currEvent = events[i];

        // first event for the cell. set the melting time and continue
        label celli = currEvent[0];

        if (currEvent[0] != prevEvent[0])
        {
            tm.setSize(pointsInCell[celli].size(), currEvent[1]);
            continue;
        }

        // extract event information
        scalar prevTime = prevEvent[1];

        scalar currTime = currEvent[1];

        List<scalar> psi0(prevEvent.size() - 2);
        List<scalar> psi(currEvent.size() - 2);

        for (int j = 0; j < psi.size(); j++)
        {
            psi0[j] = prevEvent[j + 2];
            psi[j]  = currEvent[j + 2];
        }

        int p = 0;
        for (const label& pointi : pointsInCell[celli])
        {
            scalar tp0 = Zero;
            scalar tp  = Zero;

            List<scalar> w = weights[pointi];
            
            for (int j = 0; j < psi.size(); j++)
            {
                tp0 += w[j]*psi0[j];
                tp  += w[j]*psi[j];
            }

            if ((tp <= iso_) && (tp0 > iso_))
            {
                const point& pt = positions[pointi];

                scalar m_ = min(max((iso_ - tp0)/(tp - tp0), 0), 1);

                data.append
                (
                    {
                        pt[0],
                        pt[1],
                        pt[2],
                        tm[p],
                        prevTime + m_*(currTime - prevTime),
                        (tp0 - tp) / (currTime - prevTime)
                    }
                );
            }
            else if ((tp > iso_) && (tp0 <= iso_))
            {
                scalar m_ = min(max((iso_ - tp0)/(tp - tp0), 0), 1);

                tm[p] = prevTime + m_*(currTime - prevTime);
            }

            p++;
        }
    }

    // write the events for each processor to their own file
    const fileName exacaPath
    (
        runTime_.rootPath()/runTime_.globalCaseName()/"ExaCA"
    );

    OFstream os
    (
        exacaPath + "/" + "data_" + Foam::name(Pstream::myProcNo()) + ".csv"
    );

    for(int i=0; i < data.size(); i++)
    {
        int n = data[i].size()-1;

        for(int j=0; j < n; j++)
        {
            os << data[i][j] << ",";
        }
        os << data[i][n] << "\n";
    }
}


void Foam::foamToExaCA::mapPoints(const meshSearch& searchEngine)
{
    // find event sub-space before constructing interpolants
    const pointField& points = mesh_.points();

    const vector extend = 1e-10*vector::one;

    treeBoundBox eventBb(point::max, point::min);

    for (label i = 1; i < events.size(); i++)
    {
        const label celli = events[i][0];

        if (celli == events[i - 1][0])
        {
            continue;
        }

        const labelList& vertices = mesh_.cellPoints()[celli];

        forAll(vertices, i)
        {
            eventBb.min() = min(eventBb.min(), points[vertices[i]] - extend);
            eventBb.max() = max(eventBb.max(), points[vertices[i]] + extend);
        }
    }

    label pI = 0;
    label seedi = events[0][0];

    pointsInCell.setSize(mesh_.nCells());

    const labelVector nPoints(vector::one + box_.span() / dx_);

    for (label k=0; k < nPoints.z(); ++k)
    {
        for (label j=0; j < nPoints.y(); ++j)
        {
            for (label i=0; i < nPoints.x(); ++i)
            {
                const point pt = box_.max() - vector(i, j, k)*dx_;

                if (eventBb.contains(pt))
                {
                    // shift point during search to avoid edges in pointMVC
                    const point spt = pt - vector::one*1e-10;

                    label celli = searchEngine.findCell(spt, seedi, true);

                    if (celli != -1)
                    {
                        positions.append(pt);

                        pointMVCWeight cpw(mesh_, spt, celli);

                        weights.append(cpw.weights());

                        pointsInCell[celli].append(pI);

                        pI++;
                    }
                    
                    seedi = celli;
                }
            }
        }
    }

    positions.shrink();

    weights.shrink();

    for (auto& pic : pointsInCell)
    {
        pic.shrink();
    }
}

void Foam::foamToExaCA::write()
{
    if (!execute_)
    {
        return;
    }

    // write the event data for each processor to separate files
    const fileName exacaPath
    (
        runTime_.rootPath()/runTime_.globalCaseName()/"ExaCA"
    );

    mkDir(exacaPath);

    Info<< "Number of solidification events: "
        << returnReduce(events.size(), sumOp<scalar>()) << endl;

    events.shrink();

    sort(events);

    // write exaca reduced data format with remelting
    runTime_.cpuTimeIncrement();

    meshSearch searchEngine(mesh_, polyMesh::CELL_TETS);

    if (events.size() > 0)
    {
        mapPoints(searchEngine);

        interpolatePoints();
    }

    Info<< "Successfully mapped events to exaca file in: "
        << returnReduce(runTime_.cpuTimeIncrement(), maxOp<scalar>()) << " s"
        << endl << endl;
}

// ************************************************************************* //
