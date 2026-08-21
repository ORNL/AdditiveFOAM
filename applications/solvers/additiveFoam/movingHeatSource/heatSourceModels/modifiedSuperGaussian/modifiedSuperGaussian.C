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

#include "modifiedSuperGaussian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(modifiedSuperGaussian, 0);
    addToRunTimeSelectionTable
    (
        heatSourceModel,
        modifiedSuperGaussian,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatSourceModels::modifiedSuperGaussian::modifiedSuperGaussian
(
    const word& sourceName,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(typeName, sourceName, dict, mesh),
    s_(vector::zero)
{
    k_ = heatSourceModelCoeffs_.lookup<scalar>("k");
    m_ = heatSourceModelCoeffs_.lookup<scalar>("m");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatSourceModels::modifiedSuperGaussian::update()
{
    updateDimensions();

    s_ = dimensions_/Foam::pow(2.0, 1.0/k_);
    s_.z() = dimensions_.z();

    const scalar fMax =
        Foam::pow
        (
            invIncGammaRatio_P(2.0/k_, 1.0 - profileTol_),
            1.0/k_
        );

    sourceMin_ = vector(-fMax*s_.x(), -fMax*s_.y(), -s_.z());

    sourceMax_ = vector(fMax*s_.x(), fMax*s_.y(), 0.0);

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            s_.x()*s_.y()*s_.z()*pi*Foam::tgamma(1.0 + 2.0/k_)
           *Foam::tgamma(1.0 + 1.0/m_)*Foam::tgamma(1.0 + 2.0/m_)
           /Foam::tgamma(1.0 + 3.0/m_)
        );
}


inline Foam::scalar
Foam::heatSourceModels::modifiedSuperGaussian::weight(const vector& r) const
{
    if ((r.z() > 0) || (r.z() <= -s_.z()))
    {
        return 0;
    }

    const scalar g =
        Foam::pow
        (
            1.0 - Foam::pow(-r.z()/s_.z(), m_),
            1.0/m_
        );

    const scalar f =
        Foam::pow
        (
            Foam::sqr(r.x()/(s_.x()*g))
          + Foam::sqr(r.y()/(s_.y()*g)),
            k_/2.0
        );

    return Foam::exp(-f);
}

bool Foam::heatSourceModels::modifiedSuperGaussian::read()
{
    if (heatSourceModel::read())
    {
        //- Mandatory entries
        heatSourceModelCoeffs_.lookup("k") >> k_;
        heatSourceModelCoeffs_.lookup("m") >> m_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
