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
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    h_(0.0),
    emissivity_(0.0),
    Tinf_(p.size(), Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;
}


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
    emissivity_(dict.lookup<scalar>("emissivity")),
    Tinf_("Tinf", dict, p.size())
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
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

void Foam::mixedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(Tinf_, Tinf_);
}


void Foam::mixedTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const mixedTemperatureFvPatchScalarField& tiptf =
        refCast<const mixedTemperatureFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);
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

    refValue() = Tinf_;

    const scalarField& Tp(*this);

    const scalarField& kappa_ =
        patch().lookupPatchField<volScalarField, scalar>("kappa");

    const scalar sigma_(5.67e-8);
    
    scalarField hEff_
    (
        h_ + sigma_*emissivity_*(sqr(Tp) + sqr(Tinf_))*(Tp + Tinf_)
    );

    valueFraction() = 
        1.0/
        (
            1.0
          + kappa_*patch().deltaCoeffs()/max(hEff_, 1e-15)
        );    

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
