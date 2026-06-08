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

#include "pathVector.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pathVector::pathVector()
:
    startPosition_(Zero),
    endPosition_(Zero),
    startTime_(Zero),
    endTime_(Zero),
    power_(Zero),
    displacement_(Zero),
    distance_(Zero),
    duration_(Zero),
    speed_(Zero)
{}


Foam::pathVector::pathVector
(
    const point& startPosition,
    const point& endPosition,
    const scalar startTime,
    const scalar endTime,
    const scalar power
)
:
    startPosition_(startPosition),
    endPosition_(endPosition),
    startTime_(startTime),
    endTime_(endTime),
    power_(power),
    displacement_(Zero),
    distance_(Zero),
    duration_(Zero),
    speed_(Zero)
{
    displacement_ = endPosition_ - startPosition_;

    distance_ = mag(displacement_);

    duration_ = max(endTime_ - startTime_, scalar(0));

    if (duration_ > small)
    {
        speed_ = distance_/duration_;
    }
    else
    {
        speed_ = Zero;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::point Foam::pathVector::position(const scalar time) const
{
    if (duration_ <= small)
    {
        return endPosition_;
    }

    scalar fraction = (time - startTime_)/duration_;

    fraction = min(max(fraction, scalar(0)), scalar(1));

    return startPosition_ + fraction*displacement_;
}

// ************************************************************************* //
