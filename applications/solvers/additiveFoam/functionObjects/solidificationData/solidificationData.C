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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solidificationData::solidificationData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    T_(mesh_.lookupObject<VolField<scalar>>("T")),
    R_
    (
        IOobject
        (
            "R",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::ddt(T_)
    ),
    searchEngine_(mesh_, polyMesh::CELL_TETS)
{
    read(dict);
    
    setOverlapCells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::solidificationData::~solidificationData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::solidificationData::read(const dictionary& dict)
{
    box_ = dict.lookup("box");
    isoValue_ = dict.lookup<scalar>("isoValue");
    
    return true;
}

void Foam::functionObjects::solidificationData::setOverlapCells()
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


Foam::wordList Foam::functionObjects::solidificationData::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::solidificationData::execute()
{
    //- Get current time
    const scalar& time = mesh_.time().value();

    //- Get the temperature at the previous time
    const volScalarField& T0_ = T_.oldTime();
    
    //- Calculate the thermal gradient
    const volScalarField G("G", mag(fvc::grad(T_)));
    
    forAll(overlapCells, i)
    {
        label celli = overlapCells[i];
        
        // Cooled below specified isotherm
        if ((T0_[celli] > isoValue_) && (T_[celli] <= isoValue_))
        {
            const scalar Ri = mag(R_[celli]);
            const scalar Gi = max(G[celli], small);

            const vector pt = mesh_.C()[celli];

            List<scalar> event(7);
            
            event[0] = pt[0];
            event[1] = pt[1];
            event[2] = pt[2];
            event[3] = time;
            event[4] = Ri;
            event[5] = Gi;
            event[6] = Ri / Gi;
                        
            events.append(event);
        }
    }

    //- Update cooling rate
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
    
    os << "x,y,z,t,R,G,V" << endl;

    for(int i=0; i < events.size(); i++)
    {
        int n = events[i].size() - 1;

        for(int j=0; j < n; j++)
        {
            os << events[i][j] << ",";
        }
        os << events[i][n] << "\n";
    }

    Info<< "Successfully wrote solidification data in: "
        << returnReduce(mesh_.time().cpuTimeIncrement(), maxOp<scalar>()) << " s"
        << endl << endl;

    return true;
}


bool Foam::functionObjects::solidificationData::write()
{
    return true;
}


// ************************************************************************* //
