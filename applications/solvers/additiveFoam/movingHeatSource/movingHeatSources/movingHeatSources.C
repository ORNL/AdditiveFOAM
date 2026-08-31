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

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "movingHeatSources.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingHeatSources::movingHeatSources
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (
        IOobject
        (
            "heatSourceDict",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    qDot_
    (
        IOobject
        (
            "qDot",
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPower/dimVolume, 0.0)
    ),
    refinement_(nullptr)
{
    readSources();
}


void Foam::movingHeatSources::readSources()
{
    const dictionary& sourcesDict = dict_.subDict("sources");
    const wordList names(sourcesDict.toc());
    sources_.clear();
    sources_.resize(names.size());

    forAll(sources_, i)
    {
        Info << "Adding moving heat source " << names[i] << endl;
        sources_.set
        (
            i,
            new movingHeatSource
            (
                names[i],
                sourcesDict.subDict(names[i]),
                mesh_
            )
        );
    }

    refinement_ = refinementModel::New(sources_, dict_, mesh_);
}

// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::movingHeatSources::~movingHeatSources()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::movingHeatSources::adjustDeltaT(scalar& deltaT)
{
    forAll(sources_, i)
    {
        sources_[i].beam().adjustDeltaT(deltaT);
    }
}

void Foam::movingHeatSources::update()
{
    if (dict_.readIfModified())
    {
        readSources();
    }

    qDot_ = dimensionedScalar("Zero", qDot_.dimensions(), 0.0);

    forAll(sources_, i)
    {
        if (sources_[i].beam().activePath())
        {
            sources_[i].update();

            scalar pathTime = mesh_.time().value();

            const scalar deltaT = mesh_.time().deltaTValue();

            const scalar nextTime =
                min(pathTime + deltaT, sources_[i].beam().endTime());

            const scalar beamDeltaT = sources_[i].beam().deltaT();

            volScalarField sourceQDot
            (
                IOobject
                (
                    "sourceQDot",
                    mesh_.time().name(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("Zero", qDot_.dimensions(), 0.0)
            );

            while ((nextTime - pathTime) > small)
            {
                scalar dt = min(beamDeltaT, max(0, nextTime - pathTime));

                pathTime += dt;

                sources_[i].beam().move(pathTime);

                sourceQDot += dt*sources_[i].qDot();
            }

            sourceQDot /= deltaT;

            qDot_ += sourceQDot;
        }
    }

    qDot_.correctBoundaryConditions();

    refinement_->update();
}

// ************************************************************************* //
