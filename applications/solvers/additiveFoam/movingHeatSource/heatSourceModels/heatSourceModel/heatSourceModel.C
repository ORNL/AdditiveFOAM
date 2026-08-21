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

#include "heatSourceModel.H"
#include "labelVector.H"
#include "hexMatcher.H"
#include "treeBoundBox.H"
#include "thermoPath.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatSourceModel, 0);
    defineRunTimeSelectionTable(heatSourceModel, dictionary);
}

const Foam::word Foam::heatSourceModel::heatSourceDictName
(
    "heatSourceDict"
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::heatSourceModel::createIOobject
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
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModel::heatSourceModel
(
    const word& type,
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    IOdictionary(createIOobject(dict, mesh)),

    sourceName_(sourceName),
    heatSourceDict_(dict),
    sourceDict_(heatSourceDict_.optionalSubDict(sourceName_)),
    heatSourceModelCoeffs_(sourceDict_.optionalSubDict(type + "Coeffs")),
    mesh_(mesh),
    transient_(false),
    isoValue_(great),
    D4Sigma_(0),
    minimumDepth_(0),
    depth_(0),
    isoDepth_(0),
    nPoints_(vector::one),
    tolerance_(1.0e-3),
    sourceBb_(point::zero, point::zero),
    V0_("V0", dimVolume, 0),
    absorptionModel_(nullptr),
    movingBeam_(nullptr)
{
    absorptionModel_ =
        absorptionModel::New(sourceName_, heatSourceDict_, mesh_);

    movingBeam_ =
        movingBeam::New(sourceName_, heatSourceDict_, mesh_.time());

    transient_ =
        heatSourceModelCoeffs_.lookupOrDefault<Switch>("transient", false);

    if (transient_)
    {
        isoValue_ =
            heatSourceModelCoeffs_.lookupOrDefault<scalar>
            (
                "isoValue",
                thermoPath(mesh_).liquidus()
            );
    }

    nPoints_ =
        heatSourceModelCoeffs_.lookupOrDefault<labelVector>
        (
            "nPoints",
            vector::one
        );

    tolerance_ =
        heatSourceModelCoeffs_.lookupOrDefault<scalar>
        (
            "tolerance",
            1.0e-3
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModel::updateDepth()
{
    if (!transient_)
    {
        isoDepth_ = 0;
        depth_ = minimumDepth_;
        Info << "depth: " << depth_ << endl;
        return;
    }

    const vector position = movingBeam_->position();

    // Find the maximum isotherm depth within the planar source bounds
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& cc = mesh_.C();

    scalar isoDepth = 0;

    // isocontour location evaluated linearly across faces
    for (label facei=0; facei < mesh_.nInternalFaces(); facei++)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        scalar minFace = min(T[own], T[nei]);
        scalar maxFace = max(T[own], T[nei]);

        if ((minFace < isoValue_) && (maxFace >= isoValue_))
        {
            vector d = cc[nei] - cc[own];
            const vector p =
                cc[own]
              + d*(isoValue_ - T[own])/(T[nei] - T[own])
              - position;

            if
            (
                p.x() >= sourceBb_.min().x()
             && p.x() <= sourceBb_.max().x()
             && p.y() >= sourceBb_.min().y()
             && p.y() <= sourceBb_.max().y()
             && p.z() <= 0
            )
            {
                isoDepth = max(-p.z(), isoDepth);
            }
        }
    }

    // isocontour location evaluated linearly across processor faces
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(TBf, patchi)
    {
        const fvPatchScalarField& TPf = TBf[patchi];

        const labelUList& faceCells = TPf.patch().faceCells();

        if (TPf.coupled())
        {
            const vectorField ccn
            (
                cc.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField Tn(TPf.patchNeighbourField());

            forAll(faceCells, facei)
            {
                label own = faceCells[facei];

                scalar minFace = min(T[own], Tn[facei]);
                scalar maxFace = max(T[own], Tn[facei]);

                if ((minFace < isoValue_) && (maxFace >= isoValue_))
                {
                    vector d = ccn[facei] -  cc[own];
                    const vector p =
                        cc[own]
                      + d*(isoValue_ - T[own])/(Tn[facei] - T[own])
                      - position;

                    if
                    (
                        p.x() >= sourceBb_.min().x()
                     && p.x() <= sourceBb_.max().x()
                     && p.y() >= sourceBb_.min().y()
                     && p.y() <= sourceBb_.max().y()
                     && p.z() <= 0
                    )
                    {
                        isoDepth = max(-p.z(), isoDepth);
                    }
                }
            }
        }
    }

    reduce(isoDepth, maxOp<scalar>());

    isoDepth_ = isoDepth;

    depth_ = max(minimumDepth_, isoDepth_);

    Info<< "isoDepth: " << isoDepth_
        << ", depth: " << depth_ << endl;
}

Foam::tmp<Foam::volScalarField>Foam::heatSourceModel::qDot()
{
    tmp<volScalarField> tqDot
    (
        new volScalarField
        (
            IOobject
            (
                "qDot_",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimPower/dimVolume, 0.0)
        )
    );
    volScalarField& qDot_ = tqDot.ref();

    const scalar power_ = movingBeam_->power();

    if (power_ > small)
    {
        const vector position_ = movingBeam_->position();

        const scalar a =
            (transient_ ? isoDepth_ : depth_)/(0.5*D4Sigma_);

        const dimensionedScalar absorbedPower
        (
            "etaP",
            dimPower,
            absorptionModel_->eta(a)*power_
        );

        dimensionedScalar volume = V0();

        // integrate the heat source in each overlapping cell
        volScalarField weights
        (
            IOobject
            (
                "weights",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Zero", dimless, 0.0)
        );

        treeBoundBox currentSourceBb
        (
            position_ + sourceBb_.min(),
            position_ + sourceBb_.max()
        );

        const pointField& points = mesh_.points();

        hexMatcher hex;

        forAll(mesh_.cells(), celli)
        {
            treeBoundBox cellBb(point::max, point::min);

            const labelList& vertices = mesh_.cellPoints()[celli];

            forAll(vertices, i)
            {
                cellBb.min() = min(cellBb.min(), points[vertices[i]]);
                cellBb.max() = max(cellBb.max(), points[vertices[i]]);
            }

            if (cellBb.overlaps(currentSourceBb))
            {
                if (hex.isA(mesh_, celli))
                {
                    vector dx_ =
                        cmptDivide
                        (
                            sourceBb_.span(),
                            vector(nPoints_)
                        );

                    labelVector nCellPoints =
                        max
                        (
                            cmptDivide
                            (
                                cellBb.span() + small*vector::one,
                                dx_
                            ),
                            vector::one
                        );

                    dx_ = cmptDivide(cellBb.span(), vector(nCellPoints));

                    scalar dVi = dx_.x()*dx_.y()*dx_.z();

                    scalar wi = 0.0;

                    for (label k=0; k < nCellPoints.z(); ++k)
                    {
                        for (label j=0; j < nCellPoints.y(); ++j)
                        {
                            for (label i=0; i < nCellPoints.x(); ++i)
                            {
                                const point pt
                                (
                                    cellBb.max()
                                  - cmptMultiply
                                    (
                                        vector(i + 0.5, j + 0.5, k + 0.5),
                                        dx_
                                    )
                                );

                                treeBoundBox ptBb
                                (
                                    pt - 0.5*dx_,
                                    pt + 0.5*dx_
                                );

                                // calculate weight for point in beam bound box
                                if (currentSourceBb.overlaps(ptBb))
                                {
                                    wi += weight(pt - position_)*dVi;
                                }
                            }
                        }
                    }

                    weights[celli] = wi/mesh_.V()[celli];
                }
                else
                {
                    // cell is not hexahedral, evaluate at centre
                    point d = mesh_.cellCentres()[celli] - position_;

                    weights[celli] = weight(d);
                }
            }
        }

        // stabilize numerical integration errors within 95% of applied power
        dimensionedScalar sumWeights = fvc::domainIntegrate(weights);
        scalar residual = (sumWeights/volume).value();

        if (mag(1 - residual) < 0.05)
        {
            volume = sumWeights;
        }

        qDot_ = absorbedPower * weights / volume;
    }

    return tqDot;
}


bool Foam::heatSourceModel::read()
{
    if (regIOobject::read())
    {
        sourceDict_ = optionalSubDict(sourceName_);

        heatSourceModelCoeffs_ =
            sourceDict_.optionalSubDict(type() + "Coeffs");

        transient_ =
            heatSourceModelCoeffs_.lookupOrDefault<Switch>
            (
                "transient",
                false
            );

        if (transient_)
        {
            isoValue_ =
                heatSourceModelCoeffs_.lookupOrDefault<scalar>
                (
                    "isoValue",
                    thermoPath(mesh_).liquidus()
                );
        }

        nPoints_ =
            heatSourceModelCoeffs_.lookupOrDefault<labelVector>
            (
                "nPoints",
                vector::one
            );

        tolerance_ =
            heatSourceModelCoeffs_.lookupOrDefault<scalar>
            (
                "tolerance",
                1.0e-3
            );

        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
