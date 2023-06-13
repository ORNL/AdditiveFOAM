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
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcAverage.H"
#include "DynamicList.H"

#include "meshSearch.H"
#include "labelVector.H"

#include "interpolation.H"
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

    pointsInCell(mesh_.nCells()),

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

    // read remaining fields
    box_ = this->lookup("box");
    dx_  = this->lookup<scalar>("dx");
    iso_ = this->lookup<scalar>("isotherm");

    const labelVector nPoints(vector::one + box_.span()/dx_);

    Info << "box span" << box_.span() << " nPoints: " << nPoints << endl;

    // bounding box of local processor domain
    treeBoundBox procBb(mesh_.points());

    meshSearch searchEngine(mesh_, polyMesh::CELL_TETS);

    // set point and cell addressing
    label pI = 0;
    label seedi = searchEngine.findCell(box_.min(), 0, false);

    for (label k=0; k < nPoints.z(); ++k)
    {
        for (label j=0; j < nPoints.y(); ++j)
        {
            for (label i=0; i < nPoints.x(); ++i)
            {
                const point pt = box_.max() - vector(i, j, k)*dx_;

                if (procBb.contains(pt))
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

    tm_.setSize(positions.size(), 0);

    Info<< "Assigned all points to a cell in: "
        << returnReduce(runTime_.cpuTimeIncrement(), maxOp<scalar>())
        << " s" << endl;
}

void Foam::foamToExaCA::update()
{
    if (!execute_)
    {
        return;
    }

    // interpolate temperature histories to cell vertices
    const pointScalarField Tp0_("Tp0_", Tp_);

    Tp_ = vpi_.interpolate(T_);

    // calculate the local cooling rate in each cell
    volScalarField R_
    (
        "R",
        max
        (
            -fvc::ddt(T_),
            dimensionedScalar(dimTemperature/dimTime, 0.0)
        )
    );
    
    // interpolate the average cooling rate to the cell vertices
    const pointScalarField Rp_(vpi_.interpolate(fvc::average(R_)));

    const scalar dt = runTime_.deltaTValue();
    const scalar t0 = runTime_.value() - dt;

    forAll(mesh_.cells(), cellI)
    {
        // find the cells containing the interface
        bool foundMin = false;
        bool foundMax = false;

        const labelList& vertices = mesh_.cellPoints()[cellI];

        if (pointsInCell[cellI].size())
        {
            forAll(vertices, pI)
            {
                if (Tp_[vertices[pI]] <= iso_ || Tp0_[vertices[pI]] <= iso_)
                {
                    foundMin = true;
                }

                if (Tp_[vertices[pI]] >= iso_ || Tp0_[vertices[pI]] >= iso_)
                {
                    foundMax = true;
                }

                if (foundMin && foundMax)
                {
                    break;
                }
            }
        }

        // interpolate temperature histories to exaca grid
        if (foundMin && foundMax)
        {
            for (const label& pI : pointsInCell[cellI])
            {
                const scalarField& w = weights[pI];
                
                scalar tp0 = Zero;
                scalar tp  = Zero;

                forAll(vertices, i)
                {
                    tp0 += w[i]*Tp0_[vertices[i]];
                    tp  += w[i]*Tp_[vertices[i]];
                }

                if ((tp <= iso_) && (tp0 > iso_))
                {
                    const point& pt = positions[pI];

                    scalar m_ = min(max((iso_ - tp0)/(tp - tp0), 0), 1);

                    scalar cr  = Zero;

                    forAll(vertices, i)
                    {
                        cr  += w[i]*Rp_[vertices[i]];
                    }

                    data.append
                    (
                        {
                            pt[0],
                            pt[1],
                            pt[2],
                            tm_[pI],
                            t0 + m_*dt,
                            cr
                        }
                    );
                }
                else if ((tp > iso_) && (tp0 <= iso_))
                {
                    scalar m_ = min(max((iso_ - tp0)/(tp - tp0), 0), 1);
                    tm_[pI] = t0 + m_*dt;
                }
            }
        }
    }
}

void Foam::foamToExaCA::write()
{
    if (!execute_)
    {
        return;
    }

    Info<< "Number of solidification events: "
        << returnReduce(data.size(), sumOp<scalar>()) << endl;

    const fileName exacaPath
    (
        runTime_.rootPath()/runTime_.globalCaseName()/"ExaCA"
    );

    mkDir(exacaPath);

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

// ************************************************************************* //
