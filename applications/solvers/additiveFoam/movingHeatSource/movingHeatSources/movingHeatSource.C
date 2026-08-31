/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 Oak Ridge National Laboratory
     \\/     M anipulation  |
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

#include "movingHeatSource.H"
#include "fvc.H"
#include "hexMatcher.H"
#include "treeBoundBox.H"

Foam::movingHeatSource::movingHeatSource
(
    const word& name,
    const dictionary& sourceDict,
    const fvMesh& mesh
)
:
    name_(name),
    mesh_(mesh),
    beam_(sourceDict, mesh_.time()),
    absorption_(absorptionModel::New(sourceDict.subDict("absorption"), mesh_)),
    heatSource_
    (
        heatSourceModel::New
        (
            sourceDict.subDict("heatSource"),
            mesh_
        )
    ),
    references_(sourceDict, mesh_)
{
    update();
}


void Foam::movingHeatSource::update()
{
    references_.update
    (
        heatSource_->metrics(),
        heatSource_->depth(),
        beam_.position(),
        heatSource_->bounds()
    );

    const scalar referenceDepth = references_.depth();
    const scalar sourceDepth = max(heatSource_->depth(), referenceDepth);
    heatSource_->update(sourceDepth, references_.aspectRatio());

    if (references_.isothermDepth())
    {
        Info<< "referenceDepth: " << referenceDepth
            << ", depth: " << sourceDepth << endl;
    }
}

Foam::tmp<Foam::volScalarField> Foam::movingHeatSource::qDot() const
{
    tmp<volScalarField> tqDot
    (
        new volScalarField
        (
            IOobject
            (
                "qDot_" + name_,
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPower/dimVolume, 0)
        )
    );
    volScalarField& qDot = tqDot.ref();

    if (beam_.power() <= small)
    {
        return tqDot;
    }

    const vector position = beam_.position();
    const boundBox& bounds = heatSource_->bounds();
    const labelVector& nPoints = heatSource_->nPoints();
    const dimensionedScalar absorbedPower
    (
        "etaP",
        dimPower,
        absorption_->eta(references_.aspectRatio())*beam_.power()
    );
    dimensionedScalar volume(heatSource_->V0());

    volScalarField weights
    (
        IOobject
        (
            "weights_" + name_,
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0)
    );

    const treeBoundBox currentBounds
    (
        position + bounds.min(),
        position + bounds.max()
    );
    const pointField& points = mesh_.points();
    hexMatcher hex;

    forAll(mesh_.cells(), celli)
    {
        treeBoundBox cellBounds(point::max, point::min);
        const labelList& vertices = mesh_.cellPoints()[celli];

        forAll(vertices, i)
        {
            cellBounds.min() = min(cellBounds.min(), points[vertices[i]]);
            cellBounds.max() = max(cellBounds.max(), points[vertices[i]]);
        }

        if (!cellBounds.overlaps(currentBounds))
        {
            continue;
        }

        if (!hex.isA(mesh_, celli))
        {
            weights[celli] =
                heatSource_->weight(mesh_.cellCentres()[celli] - position);
            continue;
        }

        vector dx = cmptDivide(bounds.span(), vector(nPoints));
        const labelVector nSamples
        (
            max
            (
                cmptDivide(cellBounds.span() + small*vector::one, dx),
                vector::one
            )
        );
        dx = cmptDivide(cellBounds.span(), vector(nSamples));
        const scalar dVi = dx.x()*dx.y()*dx.z();
        scalar wi = 0;

        for (label k=0; k<nSamples.z(); ++k)
        {
            for (label j=0; j<nSamples.y(); ++j)
            {
                for (label i=0; i<nSamples.x(); ++i)
                {
                    const point pt
                    (
                        cellBounds.max()
                      - cmptMultiply(vector(i + 0.5, j + 0.5, k + 0.5), dx)
                    );
                    const treeBoundBox sampleBounds
                    (
                        pt - 0.5*dx,
                        pt + 0.5*dx
                    );

                    if (currentBounds.overlaps(sampleBounds))
                    {
                        wi += heatSource_->weight(pt - position)*dVi;
                    }
                }
            }
        }

        weights[celli] = wi/mesh_.V()[celli];
    }

    const dimensionedScalar sumWeights = fvc::domainIntegrate(weights);
    const scalar residual = (sumWeights/volume).value();

    if (mag(1.0 - residual) < 0.05)
    {
        volume = sumWeights;
    }

    qDot = absorbedPower*weights/volume;
    return tqDot;
}

// ************************************************************************* //
