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
    defineTypeNameAndDebug(superGaussian, 0);
    addToRunTimeSelectionTable
    (
        heatSourceProfile,
        superGaussian,
        dictionary
    );
}
}

using Foam::constant::mathematical::pi;

Foam::scalar
Foam::heatSourceProfiles::superGaussian::coefficient
(
    const dictionary& dict,
    const scalar k
)
{
    const word definition(dict.lookupOrDefault<word>("definition", "e2"));

    if (definition == "e2")
    {
        return Foam::pow(2.0, 2.0/k);
    }
    if (definition == "secondMoment")
    {
        return
            2.0*Foam::tgamma(4.0/k)/Foam::tgamma(2.0/k);
    }

    FatalIOErrorInFunction(dict)
        << "definition must be e2 or secondMoment, found " << definition
        << exit(FatalIOError);

    return 0;
}


Foam::heatSourceProfiles::superGaussian::superGaussian
(
    const dictionary& dict,
    const fvMesh&,
    const scalar tolerance
)
:
    radius_(dict.lookup<vector2D>("radius")),
    k_(dict.lookup<scalar>("k")),
    coefficient_(coefficient(dict, k_)),
    cosTheta_(0),
    sinTheta_(0),
    metrics_(),
    bounds_(point::zero, point::zero)
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

    const scalar fMax =
        Foam::pow
        (
            invIncGammaRatio_P(2.0/k_, 1.0 - tolerance),
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

    bounds_ =
        boundBox(point(-xMax, -yMax, 0), point(xMax, yMax, 0));
}


Foam::scalar Foam::heatSourceProfiles::superGaussian::weight
(
    const scalar x,
    const scalar y
) const
{
    const scalar xr = cosTheta_*x + sinTheta_*y;
    const scalar yr = -sinTheta_*x + cosTheta_*y;
    const scalar radiusSqr =
        coefficient_
       *(sqr(xr/radius_.x()) + sqr(yr/radius_.y()));

    return Foam::exp(-Foam::pow(radiusSqr, k_/2.0));
}

// ************************************************************************* //
