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

Class
    Foam::profileMetrics

Description
    Integral and second-moment metrics of an intensity distribution in the
    beam plane. D4Sigma is ordered (major, minor).

SourceFiles
    profileMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef profileMetrics_H
#define profileMetrics_H

#include "scalar.H"
#include "symmTensor2D.H"
#include "vector2D.H"

namespace Foam
{

class profileMetrics
{
    scalar integral_;

    vector2D centroid_;

    symmTensor2D covariance_;

    vector2D D4Sigma_;

    scalar azimuth_;

public:

    profileMetrics();

    void reset
    (
        const scalar M00,
        const scalar M10,
        const scalar M01,
        const scalar M20,
        const scalar M02,
        const scalar M11
    );

    inline scalar integral() const
    {
        return integral_;
    }

    inline const vector2D& centroid() const
    {
        return centroid_;
    }

    inline const symmTensor2D& covariance() const
    {
        return covariance_;
    }

    inline const vector2D& D4Sigma() const
    {
        return D4Sigma_;
    }

    inline scalar azimuth() const
    {
        return azimuth_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
