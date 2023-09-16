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

#include "marangoniFvPatchVectorField.H"
#include "symmTransformField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchField<vector>(p, iF),
    dSigmadT_(0.0),
    Tmax_(GREAT)
{}


Foam::marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<vector>(p, iF),
    dSigmadT_(dict.lookup<scalar>("dSigmadT")),
    Tmax_(dict.lookupOrDefault<scalar>("Tmax", GREAT))
{
    evaluate();
}


Foam::marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const marangoniFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchVectorField(ptf, p, iF, mapper),
    dSigmadT_(ptf.dSigmadT_),
    Tmax_(ptf.Tmax_)
{}


Foam::marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const marangoniFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(ptf, iF),
    dSigmadT_(ptf.dSigmadT_),
    Tmax_(ptf.Tmax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::marangoniFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchVectorField::autoMap(m);
}


void Foam::marangoniFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    transformFvPatchVectorField::rmap(ptf, addr);
}


void Foam::marangoniFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    transformFvPatchVectorField::reset(ptf);
}


Foam::tmp<Foam::vectorField>
Foam::marangoniFvPatchVectorField::snGrad() const
{
    const vectorField nHat(this->patch().nf());
    vectorField pif(this->patchInternalField());

    // calculate the temperature gradient on the patch
    const volScalarField& T = 
        this->internalField().mesh().lookupObject<volScalarField>("T");

    const dimensionedScalar Tmax("Tmax", dimTemperature, Tmax_);

    const tmp<volVectorField> tgradT(fvc::grad(min(T, Tmax)));

    const vectorField tGrad(tgradT().boundaryField()[this->patch().index()]);

    // calculate the marangoni coefficient
    const dictionary& dict_ =
        db().lookupObject<IOdictionary>("transportProperties");

    dimensionedScalar mu_("mu", dimDynamicViscosity, dict_.lookup("mu"));

    scalar coeff_(dSigmadT_ / mu_.value());

    // calculate the surface normal gradient on the patch
    return
    (
        transform(I - sqr(nHat), pif) - pif
      + coeff_*transform(I - sqr(nHat), tGrad) / this->patch().deltaCoeffs()
    )*this->patch().deltaCoeffs();
}


void Foam::marangoniFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // evaluate patch field as a slip condition
    const vectorField nHat(this->patch().nf());
    vectorField pif(this->patchInternalField());

    vectorField::operator=(transform(I - sqr(nHat), pif));

    transformFvPatchField<vector>::evaluate();
}

Foam::tmp<Foam::vectorField>
Foam::marangoniFvPatchVectorField::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}

void Foam::marangoniFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchVectorField::write(os);
    writeEntry(os, "dSigmadT", dSigmadT_);
    writeEntry(os, "Tmax", Tmax_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        marangoniFvPatchVectorField
    );
}

// ************************************************************************* //
