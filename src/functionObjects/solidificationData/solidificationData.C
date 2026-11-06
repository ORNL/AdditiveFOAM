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

#include "solidificationData.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvc.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "labelVector.H"
#include "thermoPath.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(solidificationData, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidificationData,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::solidificationData::correct()
{
    setOverlapCells();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solidificationData::solidificationData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    T_(mesh_.lookupObject<volScalarField>("T")),
    R_
    (
        IOobject
        (
            "R",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::ddt(T_)
    )
{
    read(dict);

    setOverlapCells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::solidificationData::~solidificationData()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::solidificationData::read(const dictionary& dict)
{
    box_ = dict.lookup("box");

    if (!dict.readIfPresent("isoValue", isoValue_))
    {
        isoValue_ = thermoPath(mesh_).liquidus();
    }

    return true;
}


void Foam::functionObjects::solidificationData::setOverlapCells()
{
    overlapCells.clear();

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


Foam::wordList Foam::functionObjects::solidificationData::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::solidificationData::execute()
{
    const scalar& time = mesh_.time().value();

    const volScalarField& T0 = T_.oldTime();

    const volVectorField G("G", fvc::grad(T_));

    forAll(overlapCells, i)
    {
        const label celli = overlapCells[i];

        if ((T0[celli] > isoValue_) && (T_[celli] <= isoValue_))
        {
            const scalar Ri = mag(R_[celli]);
            const vector Gi = G[celli];

            const vector pt = mesh_.C()[celli];

            List<scalar> event(9);

            event[0] = pt[0];
            event[1] = pt[1];
            event[2] = pt[2];
            event[3] = time;
            event[4] = Ri;
            event[5] = Gi[0];
            event[6] = Gi[1];
            event[7] = Gi[2];
            event[8] = Ri / max(mag(Gi), small);

            events.append(event);
        }
    }

    R_ = fvc::ddt(T_);

    return true;
}


bool Foam::functionObjects::solidificationData::end()
{
    Info<< "Number of solidification events: "
        << returnReduce(events.size(), sumOp<scalar>()) << endl;

    const fileName filePath
    (
        mesh_.time().rootPath()
       /mesh_.time().globalCaseName()
       /"solidificationData"
    );

    mkDir(filePath);

    OFstream os
    (
       filePath + "/" + "data_" + Foam::name(Pstream::myProcNo()) + ".csv"
    );

    os << "x,y,z,ts,cr,gx,gy,gz,v" << endl;

    forAll(events, i)
    {
        const label n = events[i].size() - 1;

        for (label j=0; j<n; j++)
        {
            os << events[i][j] << ",";
        }

        os << events[i][n] << "\n";
    }

    Info<< "Successfully wrote solidification data in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>())
        << " s" << endl << endl;

    return true;
}


bool Foam::functionObjects::solidificationData::write()
{
    return true;
}


void Foam::functionObjects::solidificationData::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::solidificationData::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::solidificationData::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::functionObjects::solidificationData::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        correct();
    }
}


// ************************************************************************* //
