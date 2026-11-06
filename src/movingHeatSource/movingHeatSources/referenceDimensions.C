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

#include "referenceDimensions.H"
#include "PstreamReduceOps.H"
#include "thermoPath.H"
#include "volFields.H"

Foam::referenceDimensions::referenceDimensions
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    widthReference_
    (
        dict.lookupOrDefault<word>("widthReference", "D4Sigma")
    ),
    depthReference_
    (
        dict.lookupOrDefault<word>("depthReference", "constant")
    ),
    component_(areaEquivalent),
    isotherm_(0),
    width_(0),
    depth_(0)
{
    if (widthReference_ != "D4Sigma")
    {
        FatalIOErrorInFunction(dict)
            << "widthReference must be D4Sigma, found "
            << widthReference_ << exit(FatalIOError);
    }

    const word D4Sigma
    (
        dict.lookupOrDefault<word>("D4Sigma", "areaEquivalent")
    );

    if (D4Sigma == "major")
    {
        component_ = major;
    }
    else if (D4Sigma == "minor")
    {
        component_ = minor;
    }
    else if (D4Sigma != "areaEquivalent")
    {
        FatalIOErrorInFunction(dict)
            << "D4Sigma must be areaEquivalent, major, or minor, found "
            << D4Sigma << exit(FatalIOError);
    }

    if (depthReference_ == "isotherm")
    {
        isotherm_ =
            dict.lookupOrDefault<scalar>
            (
                "isotherm",
                thermoPath(mesh_).liquidus()
            );
    }
    else if (depthReference_ != "constant")
    {
        FatalIOErrorInFunction(dict)
            << "depthReference must be constant or isotherm, found "
            << depthReference_ << exit(FatalIOError);
    }
}


void Foam::referenceDimensions::update
(
    const profileMetrics& metrics,
    const scalar depth,
    const vector& position,
    const boundBox& bounds
)
{
    const vector2D& D4Sigma = metrics.D4Sigma();

    if (component_ == major)
    {
        width_ = D4Sigma.x();
    }
    else if (component_ == minor)
    {
        width_ = D4Sigma.y();
    }
    else
    {
        width_ = Foam::sqrt(D4Sigma.x()*D4Sigma.y());
    }

    if (depthReference_ == "constant")
    {
        depth_ = depth;
        return;
    }

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const volVectorField& cc = mesh_.C();
    depth_ = 0;

    for (label facei=0; facei<mesh_.nInternalFaces(); ++facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        if
        (
            min(T[own], T[nei]) < isotherm_
         && max(T[own], T[nei]) >= isotherm_
        )
        {
            const vector p =
                cc[own]
              + (cc[nei] - cc[own])
               *(isotherm_ - T[own])/(T[nei] - T[own])
              - position;

            if
            (
                p.x() >= bounds.min().x()
             && p.x() <= bounds.max().x()
             && p.y() >= bounds.min().y()
             && p.y() <= bounds.max().y()
             && p.z() <= 0
            )
            {
                depth_ = max(depth_, -p.z());
            }
        }
    }

    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(TBf, patchi)
    {
        const fvPatchScalarField& TPf = TBf[patchi];

        if (!TPf.coupled())
        {
            continue;
        }

        const labelUList& faceCells = TPf.patch().faceCells();
        const vectorField ccn
        (
            cc.boundaryField()[patchi].patchNeighbourField()
        );
        const scalarField Tn(TPf.patchNeighbourField());

        forAll(faceCells, facei)
        {
            const label own = faceCells[facei];

            if
            (
                min(T[own], Tn[facei]) < isotherm_
             && max(T[own], Tn[facei]) >= isotherm_
            )
            {
                const vector p =
                    cc[own]
                  + (ccn[facei] - cc[own])
                   *(isotherm_ - T[own])/(Tn[facei] - T[own])
                  - position;

                if
                (
                    p.x() >= bounds.min().x()
                 && p.x() <= bounds.max().x()
                 && p.y() >= bounds.min().y()
                 && p.y() <= bounds.max().y()
                 && p.z() <= 0
                )
                {
                    depth_ = max(depth_, -p.z());
                }
            }
        }
    }

    reduce(depth_, maxOp<scalar>());
}

// ************************************************************************* //
