/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "mixedTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedTemperatureFvPatchScalarField::
mixedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    h_(dict.lookup<scalar>("h")),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0.0)),
    Tinf_("Tinf", dict, p.size())
{
    refGrad() = Zero;
    valueFraction() = 0.0;

    refValue() = scalarField("Tinf", dict, p.size());
    fvPatchScalarField::operator=(refValue());
}


Foam::mixedTemperatureFvPatchScalarField::
mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    h_(ptf.h_),
    emissivity_(ptf.emissivity_),
    Tinf_(mapper(ptf.Tinf_))
{}


Foam::mixedTemperatureFvPatchScalarField::
mixedTemperatureFvPatchScalarField
(
    const mixedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    h_(ptf.h_),
    emissivity_(ptf.emissivity_),
    Tinf_(ptf.Tinf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);

    const mixedTemperatureFvPatchScalarField& tiptf =
        refCast<const mixedTemperatureFvPatchScalarField>(ptf);

    mapper(Tinf_, tiptf.Tinf_);
}


void Foam::mixedTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const mixedTemperatureFvPatchScalarField& tiptf =
        refCast<const mixedTemperatureFvPatchScalarField>(ptf);

    Tinf_.reset(tiptf.Tinf_);
}


void Foam::mixedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    mixedFvPatchScalarField::refValue() = (Tinf_);
    
    const scalarField& Tp(*this);
    
    const scalarField& kappa_ =
        patch().lookupPatchField<volScalarField, scalar>("kappa");
        
    const scalar sigma_(5.67e-8);
    
    scalarField hEff_
    (
        h_ + sigma_ * emissivity_ * (sqr(Tp) + sqr(Tinf_)) * (Tp + Tinf_)
    );
    
    valueFraction() = 
        1.0 / (1.0 + kappa_ * patch().deltaCoeffs() / max(hEff_, 1e-15));
        
    refGrad() = Zero;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::mixedTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "h", h_);
    writeEntry(os, "emissivity", emissivity_);
    writeEntry(os, "Tinf", Tinf_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mixedTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
