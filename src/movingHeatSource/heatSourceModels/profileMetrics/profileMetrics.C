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

#include "profileMetrics.H"

Foam::profileMetrics::profileMetrics()
:
    integral_(0),
    centroid_(vector2D::zero),
    covariance_(symmTensor2D::zero),
    D4Sigma_(vector2D::zero),
    azimuth_(0)
{}


void Foam::profileMetrics::reset
(
    const scalar M00,
    const scalar M10,
    const scalar M01,
    const scalar M20,
    const scalar M02,
    const scalar M11
)
{
    if (!(M00 > 0))
    {
        FatalErrorInFunction
            << "Profile integral must be positive, found " << M00
            << exit(FatalError);
    }

    integral_ = M00;

    centroid_.x() = M10/M00;

    centroid_.y() = M01/M00;

    scalar varianceX = M20/M00 - sqr(centroid_.x());

    scalar varianceY = M02/M00 - sqr(centroid_.y());

    const scalar covariance = M11/M00 - centroid_.x()*centroid_.y();

    varianceX = max(varianceX, 0.0);

    varianceY = max(varianceY, 0.0);

    covariance_ = symmTensor2D(varianceX, covariance, varianceY);

    const scalar meanVariance = 0.5*(varianceX + varianceY);

    const scalar delta =
        Foam::sqrt(0.25*sqr(varianceX - varianceY) + sqr(covariance));

    const scalar majorVariance = meanVariance + delta;

    scalar minorVariance = meanVariance - delta;

    minorVariance = max(minorVariance, 0.0);

    D4Sigma_.x() = 4.0*Foam::sqrt(majorVariance);

    D4Sigma_.y() = 4.0*Foam::sqrt(minorVariance);

    const scalar circularTolerance = rootSmall*max(meanVariance, VSMALL);

    azimuth_ =
        delta > circularTolerance
      ? 0.5*Foam::atan2(2.0*covariance, varianceX - varianceY)
      : 0.0;
}

// ************************************************************************* //
