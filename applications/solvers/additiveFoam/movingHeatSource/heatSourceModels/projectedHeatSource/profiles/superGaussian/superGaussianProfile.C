/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2026 Oak Ridge National Laboratory
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

#include "superGaussianProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace heatSourceProfiles
{
    defineTypeNameAndDebug(superGaussianProfile, 0);
    addToRunTimeSelectionTable
    (
        heatSourceProfile,
        superGaussianProfile,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;

Foam::heatSourceProfiles::superGaussianProfile::superGaussianProfile
(
    const dictionary& dict,
    const fvMesh&,
    const scalar tolerance
)
:
    dimensions_(vector2D::zero),
    k_(1),
    s_(vector2D::zero),
    integral_(0),
    D4Sigma_(0),
    profileBb_(point::zero, point::zero)
{
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");

    dimensions_ = coeffs.lookup<vector2D>("dimensions");
    k_ = coeffs.lookup<scalar>("k");

    s_ = dimensions_/Foam::pow(2.0, 1.0/k_);

    integral_ =
        pi*s_.x()*s_.y()*Foam::tgamma(1.0 + 2.0/k_);

    D4Sigma_ = D4Sigma(dimensions_.x(), dimensions_.y(), k_);

    const scalar fMax =
        Foam::pow
        (
            invIncGammaRatio_P(2.0/k_, 1.0 - tolerance),
            1.0/k_
        );

    profileBb_ =
        boundBox
        (
            point(-fMax*s_.x(), -fMax*s_.y(), 0),
            point(fMax*s_.x(), fMax*s_.y(), 0)
        );
}


Foam::scalar Foam::heatSourceProfiles::superGaussianProfile::D4Sigma
(
    const scalar x,
    const scalar y,
    const scalar k
)
{
    const scalar s = Foam::sqrt(x*y)/Foam::pow(2.0, 1.0/k);
    const scalar variance =
        0.5*Foam::tgamma(4.0/k)/Foam::tgamma(2.0/k);

    return 4.0*s*Foam::sqrt(variance);
}


Foam::scalar Foam::heatSourceProfiles::superGaussianProfile::weight
(
    const scalar x,
    const scalar y
) const
{
    return
        Foam::exp
        (
           -Foam::pow
            (
                Foam::sqr(x/s_.x()) + Foam::sqr(y/s_.y()),
                k_/2.0
            )
        );
}

// ************************************************************************* //
