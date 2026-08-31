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

#include "superGaussian.H"
#include "superGaussianProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace heatSourceModels
{
    defineTypeNameAndDebug(superGaussian, 0);
    addToRunTimeSelectionTable(heatSourceModel, superGaussian, dictionary);
}
}

using Foam::constant::mathematical::pi;

Foam::heatSourceModels::superGaussian::superGaussian
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatSourceModel(dict),
    radius_(dict.lookup<vector2D>("radius")),
    k_(dict.lookup<scalar>("k")),
    coefficient_
    (
        heatSourceProfiles::superGaussian::coefficient(dict, k_)
    ),
    cosTheta_(0),
    sinTheta_(0),
    metrics_()
{
    const scalar variance =
        0.5*Foam::tgamma(4.0/k_)/Foam::tgamma(2.0/k_);

    const scalar azimuth =
        dict.lookupOrDefault<scalar>("azimuth", 0)*pi/180.0;

    cosTheta_ = Foam::cos(azimuth);

    sinTheta_ = Foam::sin(azimuth);

    const scalar sx2 = sqr(radius_.x())/coefficient_;

    const scalar sy2 = sqr(radius_.y())/coefficient_;

    const scalar integral =
        pi*radius_.x()*radius_.y()/coefficient_
       *Foam::tgamma(1.0 + 2.0/k_);

    const scalar varianceX =
        sqr(cosTheta_)*sx2*variance
      + sqr(sinTheta_)*sy2*variance;

    const scalar varianceY =
        sqr(sinTheta_)*sx2*variance
      + sqr(cosTheta_)*sy2*variance;

    const scalar covariance =
        cosTheta_*sinTheta_*(sx2 - sy2)*variance;

    metrics_.reset
    (
        integral,
        0,
        0,
        integral*varianceX,
        integral*varianceY,
        integral*covariance
    );

    update(depth_, 0);
}


void Foam::heatSourceModels::superGaussian::update
(
    const scalar depth,
    const scalar
)
{
    sourceDepth_ = depth;

    const scalar fMax =
        Foam::pow
        (
            invIncGammaRatio_P(3.0/k_, 1.0 - tolerance_),
            1.0/k_
        );

    const scalar xMax =
        fMax/Foam::sqrt(coefficient_)
       *Foam::sqrt
        (
            sqr(radius_.x()*cosTheta_)
          + sqr(radius_.y()*sinTheta_)
        );

    const scalar yMax =
        fMax/Foam::sqrt(coefficient_)
       *Foam::sqrt
        (
            sqr(radius_.x()*sinTheta_)
          + sqr(radius_.y()*cosTheta_)
        );

    const scalar zMax = fMax*sourceDepth_/Foam::pow(2.0, 1.0/k_);

    bounds_ =
        boundBox
        (
            point(-xMax, -yMax, -zMax),
            point(xMax, yMax, 0)
        );

    V0_ =
        dimensionedScalar
        (
            "V0",
            dimVolume,
            (2.0/3.0)*pi*radius_.x()*radius_.y()*sourceDepth_
           /(coefficient_*Foam::pow(2.0, 1.0/k_))
           *Foam::tgamma(1.0 + 3.0/k_)
        );
}


Foam::scalar Foam::heatSourceModels::superGaussian::weight
(
    const vector& r
) const
{
    if (r.z() > 0)
    {
        return 0;
    }

    const scalar x = cosTheta_*r.x() + sinTheta_*r.y();

    const scalar y = -sinTheta_*r.x() + cosTheta_*r.y();

    const scalar radiusSqr =
        sqr(x/radius_.x()) + sqr(y/radius_.y());

    const scalar f =
        coefficient_*radiusSqr
      + Foam::pow(2.0, 2.0/k_)*sqr(r.z()/sourceDepth_);

    return Foam::exp(-Foam::pow(f, k_/2.0));
}

// ************************************************************************* //
