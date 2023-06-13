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

Class
    segment

Description

\*---------------------------------------------------------------------------*/

#include "segment.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Set the segment time to provided value
void Foam::segment::setTime(scalar time)
{
    time_ = time;
}

void Foam::segment::setPosition(point position)
{
    position_ = position;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// construct default segment as a zero point source
Foam::segment::segment()
:
    mode_(1),
    position_(Zero),
    power_(Zero),
    parameter_(Zero),
    time_(Zero)
{
}

// set the segement properties given a space-delimited string
Foam::segment::segment(std::string line)
{
    std::stringstream lineStream(line);
    
    lineStream
        >> mode_ 
        >> position_.x()
        >> position_.y()
        >> position_.z() 
        >> power_
        >> parameter_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
