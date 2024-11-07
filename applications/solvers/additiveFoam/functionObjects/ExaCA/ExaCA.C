/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "ExaCA.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "meshSearch.H"
#include "labelVector.H"
#include "pointMVCWeight.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

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
    T_(mesh_.lookupObject<VolField<scalar>>("T")),
    vpi_(mesh_),
    Tp_
    (
        IOobject
        (
            "Tp_",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        vpi_.interpolate(T_)
    )
{
    read(dict);
    
    setOverlapCells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ExaCA::~ExaCA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ExaCA::read(const dictionary& dict)
{
    box_ = dict.lookup("box");
    
    isoValue_ = dict.lookup<scalar>("isoValue");
    
    dx_  = dict.lookup<scalar>("dx");

    return true;
}

void Foam::functionObjects::ExaCA::setOverlapCells()
{
    // create a compact cell-stencil using the overlap sub-space
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
    return wordList::null();
}


bool Foam::functionObjects::ExaCA::execute()
{    
    const pointScalarField Tp0_("Tp0_", Tp_);

    Tp_ = vpi_.interpolate(T_);

    // capture events at interface cells: {cell id, time, vertex temperatures}
    const scalar t_ = mesh_.time().value();
    const scalar t0_ = t_ - mesh_.time().deltaTValue();
        
    forAll(overlapCells, i)
    {
        label celli = overlapCells[i];

        const labelList& vertices = mesh_.cellPoints()[celli];

        label c0 = 0;
        label c1 = 0;

        forAll(vertices, pI)
        {
            if (Tp0_[vertices[pI]] >= isoValue_)
            {
                c0++;
            }

            if (Tp_[vertices[pI]] >= isoValue_)
            {
                c1++;
            }
        }
        
        const label n = vertices.size();

        // overshoot correction: interface jumped this cell during time step
        if ( !(c0 % n) && !(c1 % n) )
        {
            if  (c0 != c1)
            {
                c0 = 0;
                c1 = 1;
            }
        }
        
        // capture solidification events
        c0 %= n;
        c1 %= n;
        
        if (c0 || c1)
        {
            List<scalar> event(n + 2);

            event[0] = celli;
        
            // add previous event
            if (c1 && !c0)
            {
                event[1] = t0_;

                forAll(vertices, pI)
                {
                    event[pI + 2] = Tp0_[vertices[pI]];
                }

                events.append(event);
            }

            // add current event
            event[1] = t_;

            forAll(vertices, pI)
            {
                event[pI + 2] = Tp_[vertices[pI]];
            }

            events.append(event);
        }
    }
    
    return true;
}

void Foam::functionObjects::ExaCA::mapPoints()
{
    mesh_.time().cpuTimeIncrement();

    if (events.size() == 0)
    {
        return;
    }
              
    // find event sub-space before constructing interpolant weights
    const pointField& points = mesh_.points();

    const vector extend = 1e-10*vector::one;

    boundBox eventBb(point::max, point::min);

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

    // find interpolant weigthts for each point
    meshSearch searchEngine(mesh_, polyMesh::CELL_TETS);  

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
    
    Info<< "Successfully mapped points to mesh in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>()) << " s"
        << endl << endl;
}

void Foam::functionObjects::ExaCA::interpolate()
{
    mesh_.time().cpuTimeIncrement();
    
    if (events.size() == 0)
    {
        return;
    }
    
    // format events in ExaCA reduced data format
    DynamicList<List<scalar>> data;

    List<scalar> tm;
    tm.setSize(pointsInCell[events[0][0]].size(), events[0][1]);

    for (label i = 1; i < events.size(); i++)
    {
        const List<scalar> prevEvent = events[i - 1];

        const List<scalar> currEvent = events[i];

        label celli = currEvent[0];

        if (currEvent[0] != prevEvent[0])
        {
            tm.setSize(pointsInCell[celli].size(), currEvent[1]);
            continue;
        }

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

            if ((tp <= isoValue_) && (tp0 > isoValue_))
            {
                const point& pt = positions[pointi];

                scalar m_ = min(max((isoValue_ - tp0)/(tp - tp0), 0), 1);

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
            else if ((tp > isoValue_) && (tp0 <= isoValue_))
            {
                scalar m_ = min(max((isoValue_ - tp0)/(tp - tp0), 0), 1);

                tm[p] = prevTime + m_*(currTime - prevTime);
            }

            p++;
        }
    }

    // write the events for each processor to their own file
    const fileName exacaPath
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()/"ExaCA"
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
    
    Info<< "Successfully interpolated and wrote ExaCA data in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>()) << " s"
        << endl << endl;
}

bool Foam::functionObjects::ExaCA::end()
{
    events.shrink();

    sort(events);

    Info<< "Number of solidification events: "
        << returnReduce(events.size(), sumOp<scalar>()) << endl;
    
    mapPoints();


    const fileName exacaPath
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()/"ExaCA"
    );

    mkDir(exacaPath);
    
    interpolate();

    return true;
}


bool Foam::functionObjects::ExaCA::write()
{
    return true;
}


// ************************************************************************* //
