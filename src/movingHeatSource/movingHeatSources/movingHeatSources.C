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
#include "fvcGrad.H"

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

template<class SourceQDot>
void Foam::movingHeatSources::integrate(const SourceQDot& sourceQDot)
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

            volScalarField qDoti
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

                qDoti += dt*sourceQDot(sources_[i]);
            }

            qDoti /= deltaT;

            qDot_ += qDoti;
        }
    }

    qDot_.correctBoundaryConditions();

    refinement_->update();
}


Foam::labelList Foam::movingHeatSources::surfaceCells
(
    const volScalarField& alpha
) const
{
    // Columns are identified by their cell-center (x, y) on a 1 nm grid
    auto column = [](const point& c)
    {
        return labelPair(round(c.x()*1e9), round(c.y()*1e9));
    };

    typedef HashTable<scalar, labelPair, labelPair::Hash<>> columnHeights;

    columnHeights zTop;
    HashTable<label, labelPair, labelPair::Hash<>> topCell;

    const volVectorField& C = mesh_.C();

    forAll(alpha, celli)
    {
        if (alpha[celli] > 0.5)
        {
            const labelPair key(column(C[celli]));
            const scalar z = C[celli].z();

            auto iter = zTop.find(key);
            if (iter == zTop.end() || z > iter())
            {
                zTop.set(key, z);
                topCell.set(key, celli);
            }
        }
    }

    columnHeights zTopGlobal(zTop);
    Pstream::mapCombineGather(zTopGlobal, maxEqOp<scalar>());
    Pstream::mapCombineScatter(zTopGlobal);

    DynamicList<label> cells(topCell.size());
    forAllConstIter(columnHeights, zTop, iter)
    {
        if (iter() == zTopGlobal[iter.key()])
        {
            cells.append(topCell[iter.key()]);
        }
    }

    return labelList(cells);
}


void Foam::movingHeatSources::update()
{
    forAll(sources_, i)
    {
        if (sources_[i].absorption().scope() == "surface")
        {
            FatalErrorInFunction
                << "Absorption model " << sources_[i].absorption().type()
                << " of source " << sources_[i].name()
                << " describes a surface hit; the source is volumetric here. "
                << "Use constant or Kelly."
                << exit(FatalError);
        }
    }

    integrate
    (
        [](const movingHeatSource& source)
        {
            return source.qDot();
        }
    );
}


void Foam::movingHeatSources::update(const volScalarField& alpha)
{
    const labelList cells(surfaceCells(alpha));

    // Incidence cosine of the vertical beam on the surface of each cell,
    // from the normal of the metal fraction
    const volVectorField gradAlpha(fvc::grad(alpha));
    scalarField cosTheta(cells.size(), 1.0);
    forAll(cells, i)
    {
        const vector& g = gradAlpha[cells[i]];
        cosTheta[i] = mag(g.z())/(mag(g) + vSmall);
    }

    forAll(sources_, i)
    {
        if (sources_[i].absorption().scope() == "volumetric")
        {
            FatalErrorInFunction
                << "Absorption model " << sources_[i].absorption().type()
                << " of source " << sources_[i].name()
                << " describes a volumetric source; the beam is absorbed at "
                << "the metal surface here. Use constant or Fresnel."
                << exit(FatalError);
        }
    }

    scalarField absorbed(sources_.size(), 0), reflected(sources_.size(), 0);
    label sourcei = 0;

    integrate
    (
        [&](const movingHeatSource& source)
        {
            // the source index: sources are integrated in order
            while (sourcei < sources_.size() && &sources_[sourcei] != &source)
            {
                sourcei++;
            }
            scalar a = 0, r = 0;
            tmp<volScalarField> tq(source.qDot(cells, cosTheta, a, r));
            if (sourcei < sources_.size())
            {
                absorbed[sourcei] += a;
                reflected[sourcei] += r;
            }
            return tq;
        }
    );

    forAll(sources_, i)
    {
        if (sources_[i].beam().activePath() && sources_[i].beam().power() > small)
        {
            // powers accumulated over the sub-steps of the path integration
            // are per call; report the mean over the step
            Info<< "Beam " << sources_[i].name() << ": absorbed "
                << absorbed[i] << " W, reflected " << reflected[i]
                << " W (sum of the path sub-steps)" << endl;
        }
    }
}

// ************************************************************************* //
